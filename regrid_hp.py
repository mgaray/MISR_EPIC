from pyresample.geometry import GridDefinition, SwathDefinition
from pyresample.kd_tree import resample_nearest
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd


def regrid_hp(da, order):
    """Regrid healpix grid to lower resolution.

    Parameters
    ----------
    da : xarray.DataArray
        Input data on healpix grid.
    order : int
        Healpix resolution level.

    Returns
    -------
    da_hp : xarray.DataArray
        Remapped data.
    """
    # Pixel numbers of source and destination hp grids
    nside = 2**order
    sratio = (da.pix.nside // nside)**2

    # Case 1: Source and destination are the same
    if sratio == 1:
        return da

    # Case 2: Destination is coarser than source, take hierarchical mean
    elif sratio > 1:
        assert da.pix.nside % nside == 0
        pix2 = np.unique(da.pix.values >> int(np.log2(sratio)))
        sub_dim = np.arange(sratio)
        da = (da.pipe(expand_pix, order)
                .pipe(xr_reshape, 'pix', pix2=pix2, sub=sub_dim)
                .mean('sub', keep_attrs=True)
                .rename(pix2='pix'))

    # Case 3: Destination is finer than source, conserve area integral of input.
    else:
        assert nside % da.pix.nside == 0
        idx = np.repeat(da.pix, sratio)
        pix2 = np.unique(da.pix.values << int(np.log2(sratio)))
        pix3 = np.stack([pix2 + i for i in range(sratio)], axis=1).ravel()
        da = da.isel(pix=idx)
        da.pix.values = pix3

    _store_metadata(da, nside, nest=True)
    return da


def infer_nside(lon, lat, order_max=12):
    """Infer healpix nside parameter whose cororderponding resolution is closes to
    input grid resolution.

    Parameters
    ----------
    lon, lat : array_like
        Arrays of longitudes and latitudes.
    order_max : int, optional
        Maximum allowed healpix order.

    Returns
    -------
    nside : int
        Inferred nside parameter.
    """
    latlon_area = np.diff(lat).mean() * np.diff(lon).mean()
    nsides = [2**r for r in range(order_max+1)]
    hp_areas = np.array([hp.nside2pixarea(n, degrees=True) for n in nsides])
    idx = np.abs(hp_areas - latlon_area).argmin()
    nside = nsides[idx]
    return nside


def xr_reshape(da, dim, **kwargs):
    """Reshape a data array.

    Examples:
    >>> dar = xr_reshape(da, 'time', year=(2017, 2018), month=np.arange(12))
    >>> dar = xr_reshape(da, 'a', b=4, c=2)
    xref: https://github.com/pydata/xarray/issues/2419

    Parameters
    ----------
    da : xarray.DataArray
        Input data.
    dim : str
        Name of dimension to squash.
    **kwargs : dict
        Key value pairs mapping new dim name to new coordinate values. Values
        can be array_like of coordinate values or a single scalar denoting size.

    Returns
    -------
    da : xarray.DataArray
        Reshaped data.
    """
    names = list(kwargs.keys())
    sizes = []
    newcoords = {}
    for k, v in kwargs.items():
        if not np.ndim(v):
            newcoords[k] = np.arange(v)
            sizes.append(int(v))
        else:
            newcoords[k] = v
            sizes.append(len(v))
    coords = {k: _do_reshape(v, dim, names, sizes) for k, v in da.coords.items()
              if k != dim}
    newcoords.update(coords)
    return _do_reshape(da, dim, names, sizes, coords=newcoords)


def _do_reshape(da, dim, names, sizes, coords=None):
    dims = da.dims
    shape = da.shape
    if dim not in dims:
        return da
    idim = dims.index(dim)
    newshape = list(shape[:idim]) + sizes + list(shape[idim+1:])
    newdims = list(dims[:idim]) + names + list(dims[idim+1:])
    data = da.data.reshape(*newshape)
    return xr.DataArray(data, dims=newdims, name=da.name, attrs=da.attrs,
                        coords=coords)


def _store_metadata(da, nside, nest=True):
    da.pix.attrs['nside'] = nside
    if nest:
        da.pix.attrs['index'] = 'NESTED'
    else:
        da.pix.attrs['index'] = 'RING'

    lon, lat = hp.pix2ang(nside, da.pix, lonlat=True, nest=nest)
    lon[lon >= 180] -= 360
    #da['lon'] = ('pix'), lon
    #da['lat'] = ('pix'), lat
    da['lon'] = lon
    da['lat'] = lat


def expand_pix(da, order):
    """Expand HEALPix pixelization to include empty pixels.

    This is the inverse operation of dropping missing values from the grid
    as it ensures enough pixel values are defined to calculate a hierarchical
    mean.

    Parameters
    ----------
    da : xarray.DataArray
        Input data on healpix grid. Should have nside and pix saved.
    order : int
        Healpix resolution level.

    Returns
    -------
    xarray.DataArray
        Expanded data.
    """
    # Make this function a no-op if input is on global grid.
    nside = 2**order
    if da.pix.size == hp.nside2npix(nside):
        return da

    # Find pixel numbers on lower order grid.
    theta, phi = hp.pix2ang(da.pix.nside, da.pix, nest=True)
    pix = np.unique(hp.ang2pix(nside, theta, phi, nest=True))
    sratio = (da.pix.nside // nside)**2
    pix2 = pix << int(np.log2(sratio))
    pix3 = np.stack([pix2 + i for i in range(sratio)], axis=1).ravel()
    return da.reindex(pix=pix3)


def latlon_to_hp_nn(da, order=None, order1=None, idx=None, pix1=None, **names):
    """Remap data from latlon grid to healpix grid via nearest-neighbors.

    Converts latlon grid to a healpix grid by subdividing the destination
    healpix grid into pixels whose areas are closest to the mean cell area
    of the latlon grid, finding assigning the nearest neighbor values from the
    latlon grid at the center of these subdivided pixels, and then averaging
    them.

    Parameters
    ----------
    da : xarray.DataArray
        Input data on cartesian grid.
    order : int, optional
        Healpix resolution level. If not given, use the order
        which has a mean pixel area that's close to the input grid.
    order1 : int, optional
        Intermmediate Healpix resolution level. If not given, use the order
        which has a mean pixel area that's close to the input grid.
    idx : array_like
        1D array with same number of elements as intermmediate healpix grid.
        Elements are indices along the flattened spatial dimension of the
        input data, which should be the nearest neighbor values.
    pix1 : array_like
        1D array with pixel numbers of intermmediate grid.

    Returns
    -------
    xarray.DataArray
        Remapped data.
    """
    # Determine nside parameters for intermmediate and output grids.
    nside1 = 2**order1 if order1 else infer_nside(da.lon, da.lat)
    nside2 = 2**order if order else nside1

    # Setup nearest neighbors
    if idx is None:
        selected, pix1 = loc_nn(da, hp.nside2order(nside1))
    else:
        selected = idx
        pix1 = np.arange(hp.nside2npix(nside1)) if pix1 is None else pix1

    # Select nearest neighbor values
    da_hp1 = da.stack(pix=da.dims).isel(pix=selected)
    da_hp1 = da_hp1.assign_coords(pix=pix1)

    # If given order is same as for intermmediate grid we are done
    if nside1 == nside2:
        _store_metadata(da_hp1, nside1, nest=True)
        return da_hp1

    return regrid_hp(da_hp1, hp.nside2order(nside2))


def loc_nn(da, order, roi=6371e3):
    """Locate the indices of the nearest neighbor points.

    Parameters
    ----------
    da : xarray.DataArray
        Input data on cartesian grid.
    order : int
        Healpix resolution level.
    roi : float, optional
        Radius of influence parameter for nearest neighbors search.
    lat_name, lon_name : str, optional
        Names of lat/lon coordinates in dataset.

    Returns
    -------
    selected : array_like
        Array containing indices of selected grid cells from input data.
    pix : array_like
        healpix pixel numbers corresponding to selected grid cells.
    """
    if np.ndim(da.lat) == 1:
        lon1, lat1 = np.meshgrid(da.lon, da.lat)
    else:
        lon1, lat1 = da.lon.values, da.lat.values
    nside = 2**order

    # Get indices of latlon points and remove missing values
    lon1, lat1 = lon1.ravel(), lat1.ravel()
    idx = np.arange(lon1.size)
    missing_mask = ~np.isnan(lon1)
    lon1, lat1 = lon1[missing_mask], lat1[missing_mask]
    idx = idx[missing_mask]
    missing_mask = ~np.isnan(lat1)
    lon1, lat1 = lon1[missing_mask], lat1[missing_mask]
    idx = idx[missing_mask]
    # Latlon values for intermediate healpix grid.
    pix = np.unique(hp.ang2pix(nside, lon1, lat1, lonlat=True, nest=True))
    lon2, lat2 = hp.pix2ang(nside, pix, lonlat=True, nest=True)
    lon2[lon2 >= 180] -= 360
    lon1[lon1 >= 180] -= 360

    # Use nearest neighbors to assign indices from lat/lon grid to pixels
    # on intermmediate hp grid.
    source_def = SwathDefinition(lon1, lat1)
    target_def = SwathDefinition(lon2, lat2)
    selected = resample_nearest(source_def, idx, target_def, roi)
    return selected, pix


def latlon_to_hp_mean(da, order=None, nest=True):
    """Remap data from latlon grid to healpix grid via simple mean.

    Converts latlon grid to a healpix grid by taking the average of all
    grid cells contained in each healpix pixel.

    Parameters
    ----------
    da : xarray.DataArray
        Input data on cartesian grid.
    order : int, optional
        Healpix resolution level. If not given, use the order
        which has a mean pixel area that's close to the input grid.
    nest : bool, optional
        If True, used NESTED (hierarchical) indexing scheme (default),
        otherwise use RING (isolatitude) indexing.

    Returns
    -------
    da_hp : xarray.DataArray
        Remapped data.
    """
    # Get healpix grid info
    nside = 2**order if order else infer_nside(da.lon, da.lat)
    lons_in, lats_in = np.meshgrid(da.lon, da.lat)
    pix_data = hp.ang2pix(nside, lons_in, lats_in, lonlat=True, nest=nest)
    pix = xr.DataArray(data=pix_data, coords=da.coords, name='pix')

    # Group all grid cells by healpix pixel number tnen take the mean
    da_hp = da.groupby(pix).mean(skipna=True).reindex(pix=np.arange(hp.nside2npix(nside)))
    _store_metadata(da_hp, nside, nest=nest)
    return da_hp


def pix_to_img(da, nx, pad=1):
    """Convert HEALPix pixelization to image array.

    Parameters
    ----------
    da : xarray.DataArray
        Input data on healpix grid. Should have nside and pix saved.
    nx : int
        Size of image in x direction.
    pad : float, optional
        Additional padding around edges of domain. Default 1 degree lat/lon.

    Returns
    -------
    xarray.DataArray
        Data converted to projection image (2D) coordinates.
    """
    latra = [da.lat.min() - pad, da.lat.max() + pad]
    lonra = [da.lon.min() - pad, da.lon.max() + pad]
    proj = hp.projector.CartesianProj(xsize=nx, flipconv='geo',
                                      latra=latra, lonra=lonra)

    # Data must be global. If not, superimpose values onto global map.
    npix = hp.nside2npix(da.pix.nside)
    if da.pix.size != npix:
        data = np.full((*da.shape[:-1], npix), np.nan)
        data[..., da.pix] = da.values
    else:
        data = da.values

    vec2pix = lambda x, y, z: hp.vec2pix(da.pix.nside, x, y, z, nest=True)
    lons, lats = proj.ij2xy()
    vec = proj.xy2vec(lons, lats)
    pix = hp.vec2pix(da.pix.nside, *vec, nest=True)
    img = data[..., pix]
    dims = list(da.dims[:-1]) + ['lat', 'lon']
    coords = dict(da.coords)
    coords.update(lat=lats[:, 0], lon=lons[0], pix=(('lat', 'lon'), pix))
    return xr.DataArray(img, dims=dims, coords=coords,
                        attrs=da.attrs, name=da.name)

def plot_hp(data, ax=None, fign=None, xsize=4096, **kwargs):
    """Visualize HEALPix map via imshow().

    Parameters
    ----------
    data : array_like
        1D array of pixel values in nested indexing order.
    ax : matplotlib.axes.Axes, optional
        Matplotlib Axes to draw onto.
    xsize : int, optional
        Size of the image, default: 4096.
    kwargs : dict, optional
        Additional keyword arguments to pass to imshow().

    Returns
    -------
    matplotlib.image.AxesImage
        HEALPix map image.
    """
    ax = plt.gca() if ax is None else ax
    m = hp.cartview(data, flip='geo', nest=True,
                    return_projected_map=True, xsize=xsize)
    m.mask = np.isnan(m)
    plt.close()
    return ax.imshow(m, extent=(-180, 180, -90, 90), **kwargs)

def calc_subgrid_residual(data1, data2):
    """Calculate subgrid residual between two healpix grids.

    Parameters
    ----------
    data1, data2 : array_like
        HEALPix grids. data1 should be higher resolution than data2 and have
        the same domain. Spatial dimension should be named 'pix'.

    Returns
    -------
    array_like
        Residual. Should have same size as data1.
    """
    # Size ratio between both grids. Should be a power of 4.
    sratio = data1.size / data2.size
    assert sratio > 1

    # Pixel number mapping from data1 to data2 via bitshift operator
    pix2 = data1.pix.values >> int(np.log2(sratio))
    residual = data1 - data2.isel(pix=pix2).assign_coords(pix=data1.pix)
    return residual

def calc_subgrid_residual_from_the_next_hp_level(data1):
    """Calculate subgrid residual of data1 after upscaling
    data1 by one HEALPix level.

    Parameters
    ----------
    data1: array_like
        HEALPix grids.Spatial dimension should be named 'pix'.

    Returns
    -------
    array_like
        Residual. Should have same size as data1.
    """
    # Pixel number mapping from data1 to the coarser level via bitshift operator
    pix2 = data1.pix.values >> 2
    data2 = xr_reshape(data1, 'pix', pix2=np.unique(pix2), sub=np.arange(4)).mean('sub')
    data2 = data2.rename(pix2='pix')
    return calc_subgrid_residual(data1, data2)

def Gorski_wavelet(residual):
    """Decompose the residual of data on HEALPix grids into
       three Gorski wavelet components.

        Parameters
        ----------
        residual: array_like
        residual values calculated on HEALPix grids. Spatial dimension should be named 'pix'.

        Returns
        -------
        array_like
        Three coefficients for Gorski wavelet at the lower hp resolutaion by 1. 
        """
    pix2 = residual.pix.values >> 2
    data_coarse = xr_reshape(residual, 'pix', pix2=np.unique(pix2), sub=np.arange(4)).mean('sub')
    data_coarse = data_coarse.rename(pix2='pix')
    residual = xr_reshape(residual, 'pix', pix2=data_coarse.pix, sub=np.arange(4))
    residual0 = residual.isel(sub=0)
    residual1 = residual.isel(sub=1)
    residual2 = residual.isel(sub=2)
    residual3 = residual.isel(sub=3)
    Gorski_1 = residual1 + residual3 - residual0 - residual2
    Gorski_2 = residual2 + residual3 - residual0 - residual1
    Gorski_3 = residual0 + residual3 - residual1 - residual2
    Gorski_power = ((Gorski_1**2+Gorski_2**2+Gorski_3**2)/3.)**0.5
    return Gorski_1, Gorski_2, Gorski_3, Gorski_power

