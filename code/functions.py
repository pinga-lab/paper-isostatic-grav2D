import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from fatiando.vis import mpl
from fatiando.gravmag import prism
from fatiando.mesher import Prism, Polygon

def bathymetry_function(A,B,c,D,yc):
    '''
    This function calculates the water bottom depth by a exponential function.
    
    Input:
    A: float - parameter related to the initial depth of the seabed.
    B: float - parameter related to the maximum depth of the seabed.
    c: float - parameter related to the slope of the seabed.
    D: float - parameter related to the horizontal location of the slope
    of the seabed.
    yc: numpy array 1D - coordinates of the center of the columns forming the
    interpretation model.
    
    Output:
    bathymetry = numpy array 1D - depth of the water bottom along the profile.
    
    '''
   
    #error messages
    
    
    #function implementation
    bathymetry = A + B*(1./(1. + np.exp(-c*(yc - D))))
    
    return bathymetry

def surface_interpolate_function(surface_picks,yc):
    '''
    This function calculates the surface depth by interpolation of picks (draw matrix).
        
    Input
    suface_picks: numpy array 2D - coordinates of the profile surface.
    yc: numpy array 1D - coordinates of the center of the columns forming the
    interpretation model.
        
    Output
    surface_interpolated: numpy array 1D - depth of the surface along the profile.
    '''
   

    #function implementation
    f = interpolate.interp1d(surface_picks[:,0], surface_picks[:,1], kind = 'linear', bounds_error = False)
    surface_interpolated = f(yc)
    
    return surface_interpolated

def moho_function(S0,dw,ds,dref,dm,dc,tw,ts):
    '''
    Calculate the relief of Moho discontinuity
    by using the Airy's isostatic model, including
    the sedimentary package.
        
    Input
	S0: numpy array 1D - compensation depth.
    dw: numpy array 1D - water density.
    ds: numpy array 1D - sediments density.
    dref: numpy array 1D - reference density.
    dm: numpy array 1D - mantle density.
	dc: numpy array 1D - crust density.
    tw: numpy array 1D - thickness of the water layer along
    the profile.
    ts: numpy array 1D - thickness of the sedimentary layer along
    the profile.
        
    Output
    S: numpy array 1D - depth of the Moho along the profile.
    '''
    
    #initialization of variables
    n = len(tw)
    shape0 = ts.shape[0]
    shape1 = ts.shape[1]
    if shape0 > shape1:
        ts = np.reshape(ts,(shape1,shape0))
        ds = np.reshape(ds,(shape1,shape0))
    
    #error messages
    assert np.min(tw) > 0.0, \
        'Bathymetry must be placed below the zero coordinate'
    assert ts.size == ds.size, \
        'ts must have the same number of elements of the ds'

    result = 0
    for (ti, di) in zip(ts, ds):
        result += (np.reshape(di,(n,1)) - dc)*np.reshape(ti,(n,1))

    aux = S0*(dref - dc)/(dm - dc)
    assert np.alltrue((tw*(dw - dc) + result)/(dm - dc) < aux), \
        'Isostatic compensation violated'

    #function implementation
    S = (tw*(dw - dc) + result + S0*(dm - dref))/(dm - dc)

    assert np.alltrue(S <= S0), \
    'S must be lower than or equal to S0'

    return S

def prism_w_function(xmax,xmin,dy,edge,dw,dref,tw,yc,border=False,last=False):
    '''
    This function create a set of 3D right rectangular prism for water layer.
        
    Input
	xmax: float - maximum coordinate of the columns forming the
    interpretation model (North direction).
    xmin: float - minimum coordinate of the columns forming the
    interpretation model (North direction).
    dy: float - Thickness of columns forming the
    interpretation model (East direction).
	edge: float - Edge extension.
    dw: numpy array 1D - water density.
    dref: numpy array 1D - reference density.
    tw: numpy array 1D - thickness of the water layer along
    the profile.
    yc: numpy array 1D - coordinates of the center of the columns forming the
    interpretation model (East direction).
	border: boolean - if True, the edge extension will be applied in the first or 
	the last prism. Otherwise, it will be applied in both.
	last: boolean - if True, the edge extension will be applied in the last prism. 
	Otherwise, it will be applied in the first.
        
    Output
    prism_w: numpy array 1D - set of 3D right rectangular prism for water layer.
    '''
    
    #initialization of variables
    n = len(yc)
    prism_w = []
    dy2 = 0.5*dy
    
    #error messages
    assert np.min(tw) > 0.0, \
        'Bathymetry must be placed below the zero coordinate'

    #function implementation

    #prisms` list of basin

    for (yi, twi) in zip(yc, tw):
    
        prism_w.append(Prism(xmin, xmax,
                             yi - dy2, yi + dy2,
                             0.0, twi, {'density': dw-dref}))
    
    if border:
        if last:
            prism_w[n-1].y2 = prism_w[n-1].y2 + edge
        else:
            prism_w[0].y1 = prism_w[0].y1 - edge
    else:
        prism_w[0].y1 = prism_w[0].y1 - edge
    	prism_w[n-1].y2 = prism_w[n-1].y2 + edge

    return prism_w

def prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,p,yc,ts0=None,ts1=None,two_layers=False,three_layers=False,border=False,last=False):
    '''
    This function create a set of 3D right rectangular prism for water layer.
        
    Input
	xmax: float - maximum coordinate of the columns forming the
    interpretation model (North direction).
    xmin: float - minimum coordinate of the columns forming the
    interpretation model (North direction).
    dy: float - Thickness of columns forming the
    interpretation model (East direction).
    edge: float - Edge extension.
    ds: numpy array 1D - sediments density.
    dref: numpy array 1D - reference density.
    tw: numpy array 1D - thickness of the water layer along
    the profile.
    p: numpy array 1D - vector parameters of the model.
    yc: numpy array 1D - coordinates of the center of the columns forming the
    interpretation model (East direction).
    ts0: numpy array 1D - if not None, thickness of the first sediment layer along
    the profile.
    ts1: numpy array 1D - if not None, thickness of the second sediment layer along
    the profile.
    two_layers: boolean - if true, the sediment layer has 2 layers. Otherwise it will has 1.
    three_layers: boolean - if true, the sediment layer has 3 layers. Otherwise it will has 1.
	border: boolean - if True, the edge extension will be applied in the first or 
	the last prism. Otherwise, it will be applied in both.
	last: boolean - if True, the edge extension will be applied in the last prism. 
	Otherwise, it will be applied in the first.

    Output
    prism_s: numpy array 1D - set of 3D right rectangular prism for sediments layer.
    '''
    
    #initialization of variables
    n = len(yc)
    prism_s = []
    prism_s0 = []
    prism_s1 = []
    prism_s2 = []
    dy2 = 0.5*dy
    
    #error messages
    assert np.min(tw) > 0.0, \
        'Bathymetry must be placed below the zero coordinate'
    
    #function implementation
    if two_layers:
        ts1 = p[0:n]
        sed0 = tw + ts0
        sed1 = tw + ts0 + ts1
        ds0 = np.copy(ds[0,:])
        ds1 = np.copy(ds[1,:])
        #prisms` list of basin
        for (yi, twi, sed0i, ds0i) in zip(yc, tw, sed0, ds0):
            prism_s0.append(Prism(xmin, xmax,
                                 yi - dy2, yi + dy2,
                                 twi, sed0i, {'density': ds0i-dref}))
        
        for (yi, sed0i, sed1i, ds1i) in zip(yc, sed0, sed1, ds1):
            prism_s1.append(Prism(xmin, xmax,
                                  yi - dy2, yi + dy2,
                                  sed0i, sed1i, {'density': ds1i-dref}))
        
        if border:
            if last:
                prism_s0[n-1].y2 = prism_s0[n-1].y2 + edge
                prism_s1[n-1].y2 = prism_s1[n-1].y2 + edge
            else:
                prism_s0[0].y1 = prism_s0[0].y1 - edge
                prism_s1[0].y1 = prism_s1[0].y1 - edge
        else:
            prism_s0[0].y1 = prism_s0[0].y1 - edge
            prism_s0[n-1].y2 = prism_s0[n-1].y2 + edge
            prism_s1[0].y1 = prism_s1[0].y1 - edge
            prism_s1[n-1].y2 = prism_s1[n-1].y2 + edge

    elif three_layers:
        ts2 = p[0:n]
        sed0 = tw + ts0
        sed1 = tw + ts0 + ts1
        sed2 = tw + ts0 + ts1 + ts2
        ds0 = np.copy(ds[0,:])
        ds1 = np.copy(ds[1,:])
        ds2 = np.copy(ds[2,:])
        #prisms` list of basin
        for (yi, twi, sed0i, ds0i) in zip(yc, tw, sed0, ds0):
            prism_s0.append(Prism(xmin, xmax,
                                  yi - dy2, yi + dy2,
                                  twi, sed0i, {'density': ds0i-dref}))
        
        for (yi, sed0i, sed1i, ds1i) in zip(yc, sed0, sed1, ds1):
            prism_s1.append(Prism(xmin, xmax,
                                  yi - dy2, yi + dy2,
                                  sed0i, sed1i, {'density': ds1i-dref}))
        
        for (yi, sed1i, sed2i, ds2i) in zip(yc, sed1, sed2, ds2):
            prism_s2.append(Prism(xmin, xmax,
                                  yi - dy2, yi + dy2,
                                  sed1i, sed2i, {'density': ds2i-dref}))
        
        if border:
            if last:
                prism_s0[n-1].y2 = prism_s0[n-1].y2 + edge
                prism_s1[n-1].y2 = prism_s1[n-1].y2 + edge
                prism_s2[n-1].y2 = prism_s2[n-1].y2 + edge
            else:
                prism_s0[0].y1 = prism_s0[0].y1 - edge
                prism_s1[0].y1 = prism_s1[0].y1 - edge
                prism_s2[0].y1 = prism_s2[0].y1 - edge
        else:
            prism_s0[0].y1 = prism_s0[0].y1 - edge
            prism_s0[n-1].y2 = prism_s0[n-1].y2 + edge
            prism_s1[0].y1 = prism_s1[0].y1 - edge
            prism_s1[n-1].y2 = prism_s1[n-1].y2 + edge
            prism_s2[0].y1 = prism_s2[0].y1 - edge
            prism_s2[n-1].y2 = prism_s2[n-1].y2 + edge

    else:
        ts = p[0:n]
        sed = tw + ts
        #prisms` list of basin
        if len(ds) < n:
            aux = ds[0]
            ds = np.empty(n)
            ds.fill(aux)
    
        for (yi, twi, sedi, dsi) in zip(yc, tw, sed, ds):
            prism_s.append(Prism(xmin, xmax,
                                 yi - dy2, yi + dy2,
                                 twi, sedi, {'density': dsi-dref}))
    
        if border:
            if last:
                prism_s[n-1].y2 = prism_s[n-1].y2 + edge
            else:
                prism_s[0].y1 = prism_s[0].y1 - edge
        else:
            prism_s[0].y1 = prism_s[0].y1 - edge
    	    prism_s[n-1].y2 = prism_s[n-1].y2 + edge
    
    return prism_s, prism_s0, prism_s1, prism_s2

def prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,p,yc,ts0=None,ts1=None,two_layers=False,three_layers=False,border=False,last=False):
    '''
    This function create a set of 3D right rectangular prism for water layer.
        
    Input
	xmax: float - maximum coordinate of the columns forming the
    interpretation model (North direction).
    xmin: float - minimum coordinate of the columns forming the
    interpretation model (North direction).
    dy: float - Thickness of columns forming the
    interpretation model (East direction).
    edge: float - Edge extension.
    S0: float - isostatic compensation surface
    dref: numpy array 1D - reference density.
	dc: numpy array 1D - crust density.
    tw: numpy array 1D - thickness of the water layer along
    the profile.
    p: numpy array 1D - vector parameters of the model.
    yc: numpy array 1D - coordinates of the center of the columns forming the
    interpretation model (East direction).
    ts0: numpy array 1D - if not None, thickness of the first sediment layer along
    the profile.
    ts1: numpy array 1D - if not None, thickness of the second sediment layer along
    the profile.
    two_layers: boolean - if true, the sediment layer has 2 layers. Otherwise it will has 1.
    three_layers: boolean - if true, the sediment layer has 3 layers. Otherwise it will has 1.
	border: boolean - if True, the edge extension will be applied in the first or 
	the last prism. Otherwise, it will be applied in both.
	last: boolean - if True, the edge extension will be applied in the last prism. 
	Otherwise, it will be applied in the first.

    Output
    prism_c: numpy array 1D - set of 3D right rectangular prism for oceanic crust layer.
    '''
    
    #initialization of variables
    n = len(yc)
    prism_c = []
    dy2 = 0.5*dy
    tm = p[n:n+n]
    dS0 = p[n+n]
    
    if two_layers:
        ts1 = p[0:n]
        sed = tw + ts0 + ts1

    elif three_layers:
        ts2 = p[0:n]
        sed = tw + ts0 + ts1 + ts2

    else:
        ts = p[0:n]
        sed = tw + ts

    #error messages
    assert np.min(tw) > 0.0, \
        'Bathymetry must be placed below the zero coordinate'

    #function implementation
    #prisms` list of basin
    for (yi, sedi, tmi, dci) in zip(yc, sed, tm, dc):
        prism_c.append(Prism(xmin, xmax,
                             yi - dy2, yi + dy2,
                             sedi, S0-tmi, {'density': dci-dref}))

    if border:
        if last:
            prism_c[n-1].y2 = prism_c[n-1].y2 + edge
        else:
            prism_c[0].y1 = prism_c[0].y1 - edge
    else:
        prism_c[0].y1 = prism_c[0].y1 - edge
    	prism_c[n-1].y2 = prism_c[n-1].y2 + edge

    return prism_c

def prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,p,yc,border=False,last=False):
    '''
    This function create a set of 3D right rectangular prism for each portion of the basin.
        
    Input
	xmax: float - maximum coordinate of the columns forming the
    interpretation model (North direction).
    xmin: float - minimum coordinate of the columns forming the
    interpretation model (North direction).
    dy: float - Thickness of columns forming the
    interpretation model (East direction).
    edge: float - Edge extension.
    S0: float - isostatic compensation surface
    dref: numpy array 1D - reference density.
    dm: numpy array 1D - mantle density.
    p: numpy array 1D - vector parameters of the model.
    yc: numpy array 1D - coordinates of the center of the columns forming the
    interpretation model (East direction).
	border: boolean - if True, the edge extension will be applied in the first or 
	the last prism. Otherwise, it will be applied in both.
	last: boolean - if True, the edge extension will be applied in the last prism. 
	Otherwise, it will be applied in the first.
        
    Output
    prism_m: numpy array 1D - set of 3D right rectangular prism for mantle layer.
    '''
    
    #initialization of variables
    n = len(yc)
    prism_m = []
    dy2 = 0.5*dy
    tm = p[n:n+n]
    dS0 = p[n+n]
    
    #error messages
   
    #function implementation
    #prisms` list of basin
    for (yi, tmi) in zip(yc, tm):
        prism_m.append(Prism(xmin, xmax,
                             yi - dy2, yi + dy2,
                             S0-tmi, S0+dS0, {'density': dm-dref}))

    if border:
        if last:
            prism_m[n-1].y2 = prism_m[n-1].y2 + edge
        else:
            prism_m[0].y1 = prism_m[0].y1 - edge
    else:
        prism_m[0].y1 = prism_m[0].y1 - edge
    	prism_m[n-1].y2 = prism_m[n-1].y2 + edge

    return prism_m

def g_function(x,yc,z,gzw,prism_s,prism_c,prism_m):
    '''
    This function calculates the gravity disturbance along the profile of the
    sedimentary basin considering a homogeneous oceanic crust.
        
    Input
    x: numpy array 1D - horizontal coordinates in the North direction along the profile.
    yc: numpy array 1D - coordinates of the center of the columns forming the
    interpretation model (East direction).
    z: numpy array 1D - vertical coordinates (Down direction) along the profile.
    prism_w: numpy array 1D - set of 3D right rectangular prism for water layer.
    prism_s: numpy array 1D - set of 3D right rectangular prism for sediments layer.
    prism_c: numpy array 1D - set of 3D right rectangular prism for oceanic crust layer.
    prism_m: numpy array 1D - set of 3D right rectangular prism for mantle layer.
 
    Output
	g: numpy array 1D - gravity disturbance along the profile.
    '''
    
    #initialization of variables
    n = len(yc)
    
    #error messages
    assert gzw.size == yc.size, \
        'gzw must have the same number of elements of the observation vector'
    
    #function implementation
    
    gzs = prism.gz(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),prism_s[0])
    gzs0 = prism.gz(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),prism_s[1])
    gzs1 = prism.gz(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),prism_s[2])
    gzs2 = prism.gz(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),prism_s[3])
    gzc = prism.gz(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),prism_c)
    gzm = prism.gz(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),prism_m)
    
    g = gzw + gzs + gzs0 + gzs1 + gzs2 + gzc + gzm
    
    return g

def density_variation(S0,dw,ds,dref,dm,tw,ts,S,middle=False):
    '''
    Calculate the oceanic crust density along the profile
    so that the Airy's isostatic model is satisfied.
        
    Input
    S0: numpy array 1D - compensation depth.
    dw: numpy array 1D - water density.
    ds: numpy array 1D - sediments density.
    dref: numpy array 1D - reference density.
    dm: numpy array 1D - mantle density.
    tw: numpy array 1D - thickness of the water layer along
    the profile.
    ts: numpy array 1D - thickness of the sedimentary layer along
    the profile.
    S: numpy array 1D - Moho depth along the profile.
    middle: boolean - if True, the density is calculated at
    the middle points. Otherwise, the density is calculated
    at the points along the profile.
        
    Output
    dc: numpy array 1D - density of the crust along the profile.
    '''

    #error messages
    assert np.max(S) <= S0, \
        'Compensation depth must be placed below the Moho'

    tc = S - tw - ts

    #function implementation
    if middle:
        dc = (S0*doc - tw*dw - ts*ds - (S0 - S)*dm)/tc

    else:
        S_profile = 0.5*(S[:-1] + S[1:])
        tw_profile = 0.5*(tw[:-1] + tw[1:])
        ts_profile = 0.5*(ts[:-1] + ts[1:])
        tc_profile = 0.5*(tc[:-1] + tc[1:])
        dc = (S0*doc - tw_profile*dw - ts_profile*ds - (S0 - S_profile)*dm)/tc_profile

    return dc

def parameter_vector(S0,prism_w,prism_s,prism_c,prism_m,two_layers=False,three_layers=False):
    '''
    This function calculates the parameter vector of the model.
        
    Input
    S0: numpy array 1D - isostatic compensation surface.
    prism_w: numpy array 1D - set of 3D right rectangular prism for water layer.
    prism_s: numpy array 1D - set of 3D right rectangular prism for sediments layer.
    prism_c: numpy array 1D - set of 3D right rectangular prism for oceanic crust layer.
    prism_m: numpy array 1D - set of 3D right rectangular prism for mantle layer.
    two_layers: boolean - if true, the sediment layer has 2 layers. Otherwise it will has 1.
    three_layers: boolean - if true, the sediment layer has 3 layers. Otherwise it will has 1.
        
    Output
    parameter_vector: numpy array 1D - parameter vector of the model.
    '''
    
    #initialization of variables
    n = len(prism_w)
    bounds_tw = []
    bounds_sed = []
    bounds_aux = []
    bounds_S = []
    bounds_SR = []
    parameter_vector = []

    #function implementation
    
    if two_layers:
        prism_aux = prism_s[1]
        prism_sed = prism_s[2]
        for i in range(n):
            bounds = prism_aux[i].get_bounds()
            bounds_aux.append(bounds[5])
        tw = np.array(bounds_aux)
    elif three_layers:
        prism_aux = prism_s[2]
        prism_sed = prism_s[3]
        for i in range(n):
            bounds = prism_aux[i].get_bounds()
            bounds_aux.append(bounds[5])
        tw = np.array(bounds_aux)
    else:
        prism_sed = prism_s[0]
        for i in range(n):
            bounds = prism_w[i].get_bounds()
            bounds_tw.append(bounds[5])
        tw = np.array(bounds_tw)

    for i in range(n):
        bounds = prism_sed[i].get_bounds()
        bounds_sed.append(bounds[5])
    sed = np.array(bounds_sed)

    ts = sed - tw

    for i in range(n):
        bounds = prism_m[i].get_bounds()
        bounds_S.append(bounds[4])
    S = np.array(bounds_S)

    bounds = prism_m[0].get_bounds()
    bounds_SR = bounds[5]
    SR = np.array(bounds_SR)

    dS0 = SR - S0
    tm = S0 - S

    parameter_vector = np.vstack((np.reshape(ts,(n,1)), np.reshape(tm,(n,1)), dS0))

    return parameter_vector

def base_known_function(dy,tw,yc,base_known,ts0=None,ts1=None,two_layers=False,three_layers=False):
    '''
    Compute known values of the parameter vector ts based on prior information.
        
    Input
    dy: float - Thickness of columns forming the interpretation model (East direction).
    tw: numpy array 1D - thickness of the water layer along the profile.
    yc: numpy array 1D - coordinates of the center of the columns forming the
    interpretation model (East direction).
    base_known: numpy array 2D - profile center coordinates and basement depths pars.
    ts0: numpy array 1D - if not None, thickness of the first sediment layer along
    the profile.
    ts1: numpy array 1D - if not None, thickness of the second sediment layer along
    the profile.
    two_layers: boolean - if true, the sediment layer has 2 layers. Otherwise it will has 1.
    three_layers: boolean - if true, the sediment layer has 3 layers. Otherwise it will has 1.
        
    Output
    r: numpy array 1D - vector of the known values parameters of the model.
    index_r: numpy array 1D - list of the index to known values parameters of the model.
        
    '''
    
    #initialization of variables
    n = len(yc)
    l = base_known.shape[0]
    r = np.zeros((l,1))
    
    base_known = base_known[np.argsort(base_known[:, 0])]
    surfaces = base_known[:,1]
    index_r = []
    indices = np.arange(n, dtype=int)
    
    if two_layers:
        saux = tw + ts0
    elif three_layers:
        saux = tw + ts0 + ts1
    else:
        saux = tw
    
    #function implementation
    
    for bi in base_known[:,0]:
        mask = np.abs(np.reshape(yc,(n,)) - bi) <= dy*0.5
        index_r.append(indices[mask])

    index_r = (np.array(index_r)).reshape((len(surfaces)))

    for i in range(l):
        r[i,:] = surfaces[i] - saux[index_r[i]]
    
    return r,index_r

def moho_known_function(dy,yc,S0,moho_known):
    '''
    Compute known values of the parameter vector based on prior information.
        
    Input
    dy: float - Thickness of columns forming the interpretation model (East direction).
    yc: numpy array 1D - coordinates of the center of the columns forming the
    interpretation model (East direction).
    S0: numpy array 1D - isostatic compensation surface.
    moho_known: numpy array 2D - profile center coordinates and moho depths pars.
        
    Output
    r: numpy array 1D - vector of the known values parameters of the model.
    index_r: numpy array 1D - list of the index to known values parameters of the model.
        
    '''
    
    #initialization of variables
    n = len(yc)
    l = moho_known.shape[0]
    r = np.zeros((l,1))
    
    moho_known = moho_known[np.argsort(moho_known[:,0])]
    surfaces = moho_known[:,1]
    index_r = []
    indices = np.arange(n, dtype=int)
    
    #function implementation
    
    for mi in moho_known[:,0]:
        mask = np.abs(np.reshape(yc,(n,)) - mi) <= dy*0.5
        index_r.append(indices[mask]+n)

    index_r = (np.array(index_r)).reshape((len(surfaces)))

    for i in range(l):
        r[i,:] = S0 - surfaces[i]
    
    return r,index_r

def T_matrix_function(pjmin,pjmax,p):
    '''
    Compute the matrix that transform the parameter vector to incorporate the inequality constraint.
        
    Input
    pjmin: float - lower limit of the interval defined for the parameter vector.
    pjmax: float - upper limit of the interval defined for the parameter vector.
    p: numpy array 1D - vector parameters of the model.
    
    Output    
    T: numpy array 2D - matrix used in the equality regularization function.
        
    '''

    #initialization of variables
    n = len(p)
    dfinv = np.zeros((n))
    
    #error messages
    #assert np.max(p) < pjmax, \
    #    'The maximum value parameter must be less than the upper limit of its interval'
    #assert np.min(p) > pjmin, \
    #    'The minimum value parameter must be greater than the lower limit of its interval'
    
    #function implementation
    dfinv = np.reshape((((pjmax - p)*(p - pjmin))/(pjmax - pjmin)),(n,))

    T = np.diag(dfinv)
    
    return T

def C_matrix_function(ds,dm,dc,two_layers=False,three_layers=False):
    '''
    Compute the diagonal matrix of mantle and unknown sediment layer density contrasts used in the Airy coinstraint function.
        
    Input
    ds: numpy array 1D - sediments density.
    dm: numpy array 1D - mantle density.
    dc: numpy array 1D - crust density.
    two_layers: boolean - if true, the sediment layer has 2 layers. Otherwise it will has 1.
    three_layers: boolean - if true, the sediment layer has 3 layers. Otherwise it will has 1.
    
    Output    
    C: numpy array 2D - diagonal matrix of mantle and unknown sediment layer density contrasts.
        
    '''
    
    #initialization of variables
    n = len(dc)
    dcaux = np.reshape(dc,(n,))
    if two_layers:
        dsaux = ds[1,:]
    elif three_layers:
        dsaux = ds[2,:]
    else:
        dsaux = np.reshape(ds,(n,))
    
    #function implementation
    cn = np.zeros((n)).reshape((n,1))
    dds = np.diag((dsaux-dcaux))
    ddm = np.diag((dm-dcaux))
    C = np.concatenate((dds, ddm, cn), axis=1)
        
    return C

def D_matrix_function(dw,dc,ds=None,two_layers=False,three_layers=False):
    '''
    Compute the diagonal matrix of water and known sediments layers density
    contrasts used in the Airy coinstraint function.
        
    Input
    dw: numpy array 1D - water density.
    dc: numpy array 1D - crust density.
    ds: numpy array 1D - sediments density.
    two_layers: boolean - if true, the sediment layer has 2 layers. Otherwise it will has 1.
    three_layers: boolean - if true, the sediment layer has 3 layers. Otherwise it will has 1.
        
    Output
    D: numpy array 2D - diagonal matrix of water and known sediments layers density
        
    '''
    
    #initialization of variables
    n = len(dc)
    dcaux = np.reshape(dc,(n,))
    
    #function implementation
    ddw = np.diag((dw-dcaux))
    if two_layers:
        ds0 = ds[0,:]
        dds0 = np.diag((ds0-dcaux))
        D = np.concatenate((ddw, dds0, dcaux.reshape((n,1))), axis=1)
    elif three_layers:
        ds0 = ds[0,:]
        ds1 = ds[1,:]
        dds0 = np.diag((ds0-dcaux))
        dds1 = np.diag((ds1-dcaux))
        D = np.concatenate((ddw, dds0, dds1, dcaux.reshape((n,1))), axis=1)
    else:
        D = np.concatenate((ddw, dcaux.reshape((n,1))), axis=1)
    
    return D

def W_matrix_function(sgm,gobs,g):
    '''
        Compute the diagonal matrix used in the Airy coinstraint function.
        
        Input
        sgm: float - exponential function decay constant.
        gobs: numpy array 1D - observed gravity disturbance along the profile.
        g: numpy array 1D - predicted gravity disturbance along the profile.
        
        Output
        W: numpy array 2D - diagonal matrix used in the Airy coinstraint function.
        
        '''
    
    #initialization of variables
    n = len(g)
    W = np.zeros((n-1,n-1))
    res_avg = np.zeros((n-1,))
    
    #function implementation
    res = gobs - g
    for i in range(n-1):
        res_avg[i] = ((res[i+1]-res[i])/2)
        W[i,i] = np.exp(-(1/sgm)*(res_avg[i]**2))
    
    return W

def R_matrix_function(n, isostatic=False):
    '''
    Compute the finite differences matrix.
        
    Input
    n: integer - used to define the matrix dimension.
        
    Output    
    R: numpy array 2D - finite differences matrix.
        
    '''
    if isostatic:
        #initialization of variables
        R = np.zeros(((n-1),n))
    
        #function implementation
        ones = np.hstack((np.ones((n-1)), 0))
        minus_ones = np.hstack((0, (-1)*np.ones((n-1))))
        P1 = np.diag(ones) + np.diag((-1)*np.ones((n-1)), k=1)
    
        R[:,:] = P1[:n-1,:]
    else:
        #initialization of variables
        R = np.zeros((2*(n-1),2*n+1))
    
        #function implementation
        ones = np.hstack((np.ones((n-1)), 0))
        minus_ones = np.hstack((0, (-1)*np.ones((n-1))))
        P1 = np.diag(ones) + np.diag((-1)*np.ones((n-1)), k=1)
        P2 = np.diag(minus_ones) + np.diag(np.ones((n-1)), k=-1)
    
        R[:n,:n] = P1
        R[n-2:2*n,n:2*n] = P2
    
    return R

def A_matrix_function(n,r,index_r):
    '''
    Compute the matrix used in the equality regularization function to basement known values.
        
    Input
    n: integer - used to define the matrix dimension.
    r: numpy array 1D - vector of the known values parameters of the model.
    index_r: numpy array 1D - list of the index to known values parameters of the model.
        
    Output
    A: numpy array 2D - matrix used in the equality regularization function.
        
    '''
    
    #initialization of variables
    m = 2*n + 1
    l = len(r)
    Aaux = np.zeros((l,m-1))
    
    #function implementation
    for (i,j) in zip (np.arange(l),index_r[:]):
        Aaux[i,j] = 1
    A = np.hstack((Aaux,np.zeros((l,1))))
    
    return A

def B_matrix_function(n,r,index_r):
    '''
    Compute the matrix used in the equality regularization function to moho known values.
        
    Input
    n: integer - used to define the matrix dimension.
    r: numpy array 1D - vector of the known values parameters of the model.
    index_r: numpy array 1D - list of the index to known values parameters of the model.
        
    Output
    B: numpy array 2D - matrix used in the equality regularization function.
        
    '''
    
    #initialization of variables
    m = 2*n + 1
    l = len(r)
    Baux = np.zeros((l,m-1))
    
    #function implementation
    for (i,j) in zip (np.arange(l),index_r[:]):
        Baux[i,j] = 1
    B = np.hstack((Baux,np.zeros((l,1))))
    
    return B

def grad_ps0_function(S0,tw,p,W,R0,C,D,ts0=None,ts1=None,two_layers=False,three_layers=False):
    '''
    Compute the gradient of the Airy constraint function.
        
    Input
    S0: numpy array 1D - isostatic compensation surface.
    tw: numpy array 1D - thickness of the water layer along the profile.
    p: numpy array 1D - vector parameters of the model.
    W: numpy array 2D - identity matrix formed by weight of data residue used in the isostatic regularization.
    R0: numpy array 2D - finite differences matrix used in the isostatic regularization.
    C: numpy array 2D - diagonal matrix of mantle and unknown sediment layer density contrasts.
    D: numpy array 2D - diagonal matrix of water and known sediments layers density.
    ts0: numpy array 1D - if not None, thickness of the first sediment layer along
    the profile.
    ts1: numpy array 1D - if not None, thickness of the second sediment layer along
    the profile.
    two_layers: boolean - if true, the sediment layer has 2 layers. Otherwise it will has 1.
    three_layers: boolean - if true, the sediment layer has 3 layers. Otherwise it will has 1.
    
    Output    
    grad_psi: numpy array 2D - gradient of the Airy constraint function.
        
    '''
    
    #initialization of variables
    n = len(tw)
    
    #error messages
    assert np.min(tw) > 0.0, \
        'Bathymetry must be placed below the zero coordinate'
    assert C.shape[0] == (p[0:n]).size, \
        'C 1st dimention must be the same of number of observation point'
    assert C.shape[1] == p.size, \
        'C 2nd dimention must be the same of number of parameter vector elements'

    if two_layers:
        t = np.vstack((tw, ts0, S0))
    elif three_layers:
        t = np.vstack((tw, ts0, ts1, S0))
    else:
        t = np.vstack((tw, S0))
    
    #function implementation
    grad_psi = 2.*(((((C.T.dot(R0.T)).dot(W.T)).dot(W)).dot(R0)).dot((D.dot(t) + C.dot(p))))
    
    return grad_psi

def grad_psi1_function(p,R):
    '''
    Compute gradient of the first order Tikhonov regularization function.
        
    Input
    p: numpy array 1D - vector parameters of the model.
    R: numpy array 2D - finite differences matrix used in the first order Tikhonov regularization.
        
    Output    
    grad_psi: numpy array 1D - gradient of the first order Tikhonov regularization function.
        
    '''
    
    #error messages
    assert R.shape[0] == 2*(int((len(p)-1)*0.5)-1), \
       'R 1st dimention must be the same that 2*(n-1)'
    assert R.shape[1] == p.size, \
       'R 2nd dimention must be the same of number of parameter vector elements'
    
    
    #function implementation
    grad_psi = 2.*(R.T.dot(R)).dot(p)
    
    return grad_psi

def grad_psi2_function(p,r,Matrix):
    '''
    Compute gradient of the equality constraint function.
        
    Input
    p: numpy array 1D - vector parameters of the model.
    r: numpy array 1D - vector of the known values parameters of the model.
    Matrix: numpy array 2D - matrix used in the equality regularization function.
        
    Output
    grad_psi: numpy array 1D - gradient of the equality constraint function.
        
    '''
    
    #error messages
    assert Matrix.shape[0] == r.size, \
        'Matrix 1nd dimention must be the same of number of moho known values vector elements'
    assert Matrix.shape[1] == p.size, \
        'Matrix 2nd dimention must be the same of number of parameter vector elements'

    #function implementation
    grad_psi = 2.*Matrix.T.dot(Matrix.dot(p)-r)

    return grad_psi

def gama_function(alpha0,alpha1,alpha2,alpha3,lamb,S0,tw,gobs,g,p,rs,rm,W,R0,C,D,R,A,B,ts0=None,ts1=None,two_layers=False,three_layers=False):
    
    '''
    Compute sum of residue and constraint functions.
        
    Input
    alpha0: float - weight of parameter of regularization (isostatic constrain).
    alpha1: float - weight of parameter of regularization (TK1 constrain).
    alpha2: float - weight of parameter of regularization (equality constrain).
    alpha3: float - weight of parameter of regularization (equality constrain).
    lamb: float - parameter of regularization
    S0: numpy array 1D - isostatic compensation surface.
    tw: numpy array 1D - thickness of the water layer along the profile.
    gobs: numpy array 1D - observed gravity disturbance along the profile.
    g: numpy array 1D - predicted gravity disturbance along the profile.
    p: numpy array 1D - vector parameters of the model.
    rs: numpy array 1D - vector of the known values parameters of the model (top of basement known depths).
    rm: numpy array 1D - vector of the known values parameters of the model (Moho known depths).
    W: numpy array 2D - identity matrix formed by weight of data residue used in the isostatic regularization.
    R0: numpy array 2D - finite differences matrix used in the isostatic regularization.
    C: numpy array 2D - diagonal matrix of mantle and unknown sediment layer density contrasts.
    D: numpy array 2D - diagonal matrix of water and known sediments layers density.
    R: numpy array 2D - finite differences matrix used in the first order Tikhonov regularization.
    A: numpy array 2D - matrix used in the equality regularization function (top of basement known depths).
    B: numpy array 2D - matrix used in the equality regularization function (Moho known depths).
    ts0: numpy array 1D - if not None, thickness of the first sediment layer along
    the profile.
    ts1: numpy array 1D - if not None, thickness of the second sediment layer along
    the profile.
    two_layers: boolean - if true, the sediment layer has 2 layers. Otherwise it will has 1.
    three_layers: boolean - if true, the sediment layer has 3 layers. Otherwise it will has 1.
        
    Output
    gama: numpy array 1D - sum of residue and constraint functions.
        
    '''
    
    #initialization of variables
    n = len(gobs)
    
    #error messages
    assert gobs.size == g.size == (p[0:n]).size, \
        ' predicted data 1 and predicted data 2 must have the same number of elements of the observation vector'
    assert C.shape[0] == (p[0:n]).size, \
        'C 1st dimention must be the same of number of observation point'
    assert C.shape[1] == p.size, \
        'C 2nd dimention must be the same of number of parameter vector elements'
    assert A.shape[0] == rs.size, \
        'A 1nd dimention must be the same of number of basement known values vector elements'
    assert A.shape[1] == p.size, \
        'A 2nd dimentions must be the same of number of parameter vector elements'
    assert B.shape[0] == rm.size, \
        'B 1nd dimention must be the same of number of moho known values vector elements'
    assert B.shape[1] == p.size, \
        'B 2nd dimention must be the same of number of parameter vector elements'
    assert R.shape[0] == 2*(int((len(p)-1)*0.5)-1), \
        'R 1st dimention must be the same that 2*(n-1)'
    assert R.shape[1] == p.size, \
        'R 2nd dimention must be the same of number of parameter vector elements'
    
    #function implementation
    if two_layers:
        t = np.vstack((tw, ts0, S0))
    elif three_layers:
        t = np.vstack((tw, ts0, ts1, S0))
    else:
        t = np.vstack((tw, S0))
    
    phi = (1./n)*((gobs - g).T.dot(gobs - g))[0,0]
    psi0 = ((W.dot(R0.dot(D.dot(t))) + W.dot(R0.dot(C.dot(p)))).T.dot(W.dot(R0.dot(D.dot(t))) + W.dot(R0.dot(C.dot(p)))))[0,0]
    psi1 = ((R.dot(p)).T.dot(R.dot(p)))[0,0]
    psi2 = ((A.dot(p) - rs).T.dot((A.dot(p) - rs)))[0,0]
    psi3 = ((B.dot(p) - rm).T.dot((B.dot(p) - rm)))[0,0]

    gama = phi + lamb*(alpha0*psi0 + alpha1*psi1 + alpha2*psi2 + alpha3*psi3)

    return gama

def convergence_function(gama0,gama,beta):
    '''
    Test convergence between objective functions.
        
    Input
    gama0: float - objective function in previous iteration.
    gama: float - objective function in current iteration.
    beta: float - optimization variable that controls the convergence of the inversion.
        
    Output
    tau: Boolean - true or false.
        
    '''
    
    #function implementation
    tau = abs((gama0-gama))/gama0 < beta
    
    return tau   

def polygons_display(surfu,surfd,dS0u,dS0d,yli,yri,dy,S0,dw,ds,dm,dref,dc,tw,p,yc,area,ts0=None,ts1=None,two_layers=False,three_layers=False):
    '''
    This function display the models of prisms with a parameter and its finite difference step.
        
    Input
    surfu: float - surface chosen for viewing with subtraction of FD step.
    surfd: float - surface chosen for viewing with addition of FD step.
    dS0u: float -  difference between isostatic compensation surface and reference Moho
    with subtraction of FD step, if necessary.
    dS0d: float -  difference between isostatic compensation surface and reference Moho
    with addition of FD step, if necessary.
    yli:  float - left coordinate of the prism column relative to the chosen surface.
    yri: float - right coordinate of the prism column relative to the chosen surface.
    dy: float - Thickness of columns forming the
    S0: numpy array 1D - isostatic compensation surface.
    interpretation model (East direction).
    dw: numpy array 1D - water density.
    ds: numpy array 1D - sediments density.
    dm: numpy array 1D - mantle density.
    dref: numpy array 1D - reference density.
    dc: numpy array 1D - crust density.
    tw: numpy array 1D - thickness of the water layer along
    the profile.
    p: numpy array 1D - vector parameters of the model.
    yc: numpy array 1D - coordinates of the center of the columns forming the
    interpretation model (East direction).
    area: numpy array 2D - model  limits = [ymin, ymax, zmax, zmin]
    ts0: numpy array 1D - if not None, thickness of the first sediment layer along
    the profile.
    ts1: numpy array 1D - if not None, thickness of the second sediment layer along
    the profile.
    two_layers: boolean - if true, the sediment layer has 2 layers. Otherwise it will has 1.
    three_layers: boolean - if true, the sediment layer has 3 layers. Otherwise it will has 1.
       
    Output
    image of the prism models.
    '''
    
    #initialization of variables
    ymin = area[0]
    ymax = area[1]
    zmin = area[3]
    zmax = area[2]
    ycmin = ymin + 0.5*dy
    ycmax = ymax - 0.5*dy
    n = len(yc)
    ts = p[0:n]
    tm = p[n:n+n]
    dS0 = p[n+n]
    dy2 = 0.5*dy
    S = S0 - tm
    if two_layers:
        ts1 = np.copy(ts)
        ds0 = ds[0,:1]
        ds1 = ds[1,:1]
        ds2 = np.zeros((1))
    elif three_layers:
        ts2 = np.copy(ts)
        ds0 = ds[0,:1]
        ds1 = ds[1,:1]
        ds2 = ds[2,:1]
    else:
        ts0 = np.copy(ts)
        ds0 = ds[:1,0]
        ds1 = np.zeros((1))
        ds2 = np.zeros((1))
    
    #error messages
    assert surfu < surfd, \
        'surface plus FD step must be placed below surface less FD step'
    
    #function implementation
    polygons_water = []
    for (yi, twi) in zip(yc, tw):
        y1 = yi - dy2
        y2 = yi + dy2
        polygons_water.append(Polygon(np.array([[y1, y2, y2, y1], 
                                                [0.0, 0.0, twi, twi]]).T,
                                                props={'density': dw - dref}))
    polygons_sed0 = []
    sed0 = tw + ts0
    for (yi, twi, s0i) in zip(yc, np.reshape(tw,(n,)), np.reshape(sed0,(n,))):
        y1 = yi - dy2
        y2 = yi + dy2
        polygons_sed0.append(Polygon(np.array([[y1, y2, y2, y1],
                                              [twi, twi, s0i, s0i]]).T,
                                              props={'density': ds0 - dref}))
    if two_layers:
        polygons_sed1 = []
        polygons_sed2 = []
        sed1 = tw + ts0 + ts1
        for (yi, s0i, s1i) in zip(yc, np.reshape(sed0,(n,)), np.reshape(sed1,(n,))):
            y1 = yi - 0.5*dy
            y2 = yi + 0.5*dy
            polygons_sed1.append(Polygon(np.array([[y1, y2, y2, y1],
                                                   [s0i, s0i, s1i, s1i]]).T,
                                                   props={'density': ds1 - dref}))
        sed = np.copy(sed1)
    elif three_layers:
        polygons_sed1 = []
        sed1 = tw + ts0 + ts1
        for (yi, s0i, s1i) in zip(yc, np.reshape(sed0,(n,)), np.reshape(sed1,(n,))):
            y1 = yi - 0.5*dy
            y2 = yi + 0.5*dy
            polygons_sed1.append(Polygon(np.array([[y1, y2, y2, y1],
                                                   [s0i, s0i, s1i, s1i]]).T,
                                                   props={'density': ds1 - dref}))
        polygons_sed2 = []
        sed2 = tw + ts0 + ts1 + ts2
        for (yi, s1i, s2i) in zip(yc, np.reshape(sed1,(n,)), np.reshape(sed2,(n,))):
            y1 = yi - 0.5*dy
            y2 = yi + 0.5*dy
            polygons_sed2.append(Polygon(np.array([[y1, y2, y2, y1],
                                                   [s1i, s1i, s2i, s2i]]).T,
                                                   props={'density': ds2 - dref}))
        sed = np.copy(sed2)
    else:
        sed = np.copy(sed0)
        polygons_sed1 = []
        polygons_sed2 = []

    polygons_crust = []
    for (yi, si, Si, dci) in zip(yc, np.reshape(sed,(n,)), np.reshape(S,(n,)), dc):
        y1 = yi - dy2
        y2 = yi + dy2
        polygons_crust.append(Polygon(np.array([[y1, y2, y2, y1],
                                                [si, si, Si, Si]]).T,
                                      props={'density': dci - dref}))
    polygons_mantle = []
    for (yi, Si) in zip(yc, np.reshape(S,(n,))):
        y1 = yi - dy2
        y2 = yi + dy2
        polygons_mantle.append(Polygon(np.array([[y1, y2, y2, y1], 
                                                 [Si, Si, S0+2*dS0, S0+2*dS0]]).T,
                                                 props={'density': dm - dref}))
    #plot
    mpl.close('all')
    mpl.figure(figsize=(15,10))

    mpl.title('Sedimentary basin', fontsize=18)

    colors = np.vstack((dw, ds0, ds1, ds2, dm, dc))
    levels = plt.contourf(np.reshape(colors, (colors.size, 1)), 20, 
                          cmap=mpl.get_cmap('Greys'))
    cb = mpl.colorbar(orientation='horizontal', format='%.0f', pad=0.2)
    cb.set_label('Density (kg/m$^{3}$)', fontsize=16)
    dmax = np.max(levels.levels)
    dmin = np.min(levels.levels)

    cw = np.empty(n)
    cw.fill(((dmax - dw[0])/(dmax - dmin)))
    cs0 = np.empty(n)
    cs0.fill(((dmax - ds0[0])/(dmax - dmin)))
    cs1 = np.empty(n)
    cs1.fill(((dmax - ds1[0])/(dmax - dmin)))
    cs2 = np.empty(n)
    cs2.fill(((dmax - ds2[0])/(dmax - dmin)))
    cc = ((dmax - dc)/(dmax - dmin))
    cm = np.empty(n)
    cm.fill(((dmax - dm[0])/(dmax - dmin)))

    mpl.paths([[ymin, 0.0]], [[ymax, 0.0]], style='-k', linewidth=1)
    for (pwi, cwi) in zip(polygons_water, cw):
        mpl.polygon(pwi, style='-w', linewidth=0, fill=(cwi, cwi, cwi))
    for (ps0i, cs0i) in zip(polygons_sed0, cs0):
        mpl.polygon(ps0i, style='-w', linewidth=0, fill=(cs0i, cs0i, cs0i))
    for (ps1i, cs1i) in zip(polygons_sed1, cs1):
        mpl.polygon(ps1i, style='-w', linewidth=0, fill=(cs1i, cs1i, cs1i))
    for (ps2i, cs2i) in zip(polygons_sed2, cs2):
        mpl.polygon(ps2i, style='-w', linewidth=0, fill=(cs2i, cs2i, cs2i))
    for (pci, cci) in zip(polygons_crust, np.reshape(cc,(n,))):
        mpl.polygon(pci, style='-w', linewidth=0, fill=(cci, cci, cci))
    for (pmi, cmi) in zip(polygons_mantle, cm):
        mpl.polygon(pmi, style='-w', linewidth=0, fill=(cmi, cmi, cmi))
            
    mpl.paths([[yli, surfu]], [[yri, surfu]], style='--r', linewidth=1, label='p + dp')
    mpl.paths([[yli, surfd]], [[yri, surfd]], style='--b', linewidth=1, label='p - dp')
    #mpl.paths([[yli, S0+dS0u]], [[yri, S0+dS0u]], style='--r', linewidth=1)
    #mpl.paths([[yli, S0+dS0d]], [[yri, S0+dS0d]], style='--b', linewidth=1)
    mpl.paths([[ycmin, S0]], [[ycmax, S0]], style='--y', linewidth=1, label='isostatic compensation surface')
    mpl.paths([[ycmin, S0+dS0]], [[ycmax, S0+dS0]], style='--w', linewidth=1, label='Moho reference surface')
    mpl.ylim(S0+2*dS0, zmin)
    mpl.xlim(ycmin, ycmax)
    mpl.xlabel('y (km)', fontsize=16)
    mpl.ylabel('z (km)', fontsize=16)
    mpl.xticks(fontsize=14)
    mpl.yticks(fontsize=14)
    mpl.m2km()
    mpl.tight_layout()
    mpl.legend(loc='lower left', fontsize='small')
    mpl.show() 

    return

def G_matrix_function(xmax,xmin,dy,edge,dp1,dp2,S0,dw,ds,dm,dref,dc,tw,p,yc,ts0=None,ts1=None,two_layers=False,three_layers=False):

    '''
    This function create the sensibility matrix considering only the contribution
    of the prism column containing the disturbed parameter.
        
    Input
    xmax: float - maximum coordinate of the columns forming the
    interpretation model (North direction).
    xmin: float - minimum coordinate of the columns forming the
    interpretation model (North direction).
    dy: float - Thickness of columns forming the
    interpretation model (East direction).
    edge: float - Edge extension.
    dp1: float - ts and tm parameter variation to Finite Difference
    dp2: float - S0 and dS0 parameter variation to Finite Difference
    S0: numpy array 1D - isostatic compensation surface.
    dw: float - water density.
    ds: float - sediments density.
    dm: float - mantle density.
    dref: numpy array 1D - reference density.
    dc: numpy array 1D - crust density.
    tw: numpy array 1D - thickness of the water layer along
    the profile.
    p: numpy array 1D - vector parameters of the model.
    yc: numpy array 1D - coordinates of the center of the columns forming the
    interpretation model (East direction).
    ts0: numpy array 1D - if not None, thickness of the first sediment layer along
    the profile.
    ts1: numpy array 1D - if not None, thickness of the second sediment layer along
    the profile.
    two_layers: boolean - if true, the sediment layer has 2 layers. Otherwise it will has 1.
    three_layers: boolean - if true, the sediment layer has 3 layers. Otherwise it will has 1.
        
    Output
    G: numpy array 2D - sensibility matriz.
    '''
    
    #initialization of variables
    x = np.zeros_like(yc)
    z = np.zeros_like(yc)-150.0
    n = len(yc)
    
    dS0 = p[n+n]
    
    gplus = np.zeros_like(yc)
    gminus = np.zeros_like(yc)
    pplus = np.copy(p)
    pminus = np.copy(p)

    G = np.zeros((n,2*n+1))
    
    
    #error messages
    assert np.min(tw) > 0.0, \
        'Bathymetry must be placed below the zero coordinate'
    
    #function implementation
    prism_w = prism_w_function(xmax,xmin,dy,edge,dw,dref,tw,yc)
    gzw = prism.gz(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),prism_w)

    if two_layers:
        tc = S0 - tw - ts0 - p[0:n] - p[n:n+n]
    elif three_layers:
        tc = S0 - tw - ts0 - ts1 - p[0:n] - p[n:n+n]
    else:
        tc = S0 - tw - p[0:n] - p[n:n+n]
    
    save_edge = edge
    
    for j in range(2*n):
        dp = dp1
        if tc[j%n] <= dp:
            if j < n:
                pplus[j] = p[j]
                pminus[j] = p[j]-dp
            else:
                pplus[j] = p[j]+dp
                pminus[j] = p[j]   
        
        elif p[j] <= dp:
            if j < n:
                pplus[j] = p[j]+dp
                pminus[j] = p[j]
            else:
                pplus[j] = p[j]
                pminus[j] = p[j]-dp   
       
        else:
            pplus[j] = p[j]+dp
            pminus[j] = p[j]-dp 
        
        ts_plus = pplus[0:n]
        tm_plus = pplus[n:n+n]
        ts_minus = pminus[0:n]
        tm_minus = pminus[n:n+n]
        
        k = j%n
        
        yck = np.array([yc[k]])
        twk = np.array([tw[k]])
        dck = np.array([dc[k]])
        ts_plusk = np.array([ts_plus[k]])
        tm_plusk = np.array([tm_plus[k]])
        ts_minusk = np.array([ts_minus[k]])
        tm_minusk = np.array([tm_minus[k]])
       
        pplusk = np.vstack((ts_plusk, tm_plusk, dS0))
        pminusk = np.vstack((ts_minusk, tm_minusk, dS0))
        
        if k == 0:
            if two_layers:
                prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pplusk,yck,ts0,two_layers=True,border=True)
                prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pplusk,yck,ts0,two_layers=True,border=True)
                prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pminusk,yck,ts0,two_layers=True,border=True)
                prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pminusk,yck,ts0,two_layers=True,border=True)
            elif three_layers:
                prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pplusk,yck,ts0,ts1,three_layers=True,border=True)
                prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pplusk,yck,ts0,ts1,three_layers=True,border=True)
                prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pminusk,yck,ts0,ts1,three_layers=True,border=True)
                prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pminusk,yck,ts0,ts1,three_layers=True,border=True)
            else:
                prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pplusk,yck,border=True)
                prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pplusk,yck,border=True)
                prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pminusk,yck,border=True)
                prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pminusk,yck,border=True)

            prism_m_plus = prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,pplusk,yck,border=True)
            prism_m_minus = prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,pminusk, yck,border=True)

        elif k == (n-1):
            if two_layers:
                prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pplusk,yck,ts0,two_layers=True,border=True,last=True)
                prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pplusk,yck,ts0,two_layers=True,border=True,last=True)
                prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pminusk,yck,ts0,two_layers=True,border=True,last=True)
                prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pminusk,yck,ts0,two_layers=True,border=True,last=True)
            elif three_layers:
                prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pplusk,yck,ts0,ts1,three_layers=True,border=True,last=True)
                prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pplusk,yck,ts0,ts1,three_layers=True,border=True,last=True)
                prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pminusk,yck,ts0,ts1,three_layers=True,border=True,last=True)
                prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pminusk,yck,ts0,ts1,three_layers=True,border=True,last=True)
            else:
                prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pplusk,yck,border=True,last=True)
                prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pplusk,yck,border=True,last=True)
                prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pminusk,yck,border=True,last=True)
                prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pminusk,yck,border=True,last=True)
            
            prism_m_plus = prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,pplusk,yck,border=True,last=True)
            prism_m_minus = prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,pminusk, yck,border=True,last=True)

        else:
            edge=0.0
            if two_layers:
                prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pplusk,yck,ts0,two_layers=True)
                prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pplusk,yck,ts0,two_layers=True)
                prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pminusk,yck,ts0,two_layers=True)
                prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pminusk,yck,ts0,two_layers=True)
            elif three_layers:
                prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pplusk,yck,ts0,ts1,three_layers=True)
                prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pplusk,yck,ts0,ts1,three_layers=True)
                prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pminusk,yck,ts0,ts1,three_layers=True)
                prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pminusk,yck,ts0,ts1,three_layers=True)
            else:
                prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pplusk,yck)
                prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pplusk,yck)
                prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,twk,pminusk,yck)
                prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dck,twk,pminusk,yck)
            
            prism_m_plus = prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,pplusk,yck)
            prism_m_minus = prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,pminusk,yck)

        edge = save_edge
        
        pplus = np.copy(p)
        pminus = np.copy(p)
    
        for i in range(n):
        
            xi = np.array([x[i]])
            yci = np.array([yc[i]])
            zi = np.array([z[i]])
            gzwi = np.array([gzw[i]])
            
            gplus[i] = g_function(xi,yci,zi,gzwi,prism_s_plus,prism_c_plus,prism_m_plus)
            
            gminus[i] = g_function(xi,yci,zi,gzwi,prism_s_minus,prism_c_minus,prism_m_minus)
            
            G[i,j] = (gplus[i] - gminus[i])/(2*dp)
    
    # Para os parametros S0 e dS0:
    for k in range(2*n,2*n+1):
        dp = dp2
        pplus[k] = p[k] + dp
        pminus[k] = p[k] - dp
                  
        if two_layers:
            prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pplus,yc,ts0,two_layers=True)
            prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pplus,yc,ts0,two_layers=True)
            prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pminus,yc,ts0,two_layers=True)
            prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pminus,yc,ts0,two_layers=True)
        elif three_layers:
            prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pplus,yc,ts0,ts1,three_layers=True)
            prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pplus,yc,ts0,ts1,three_layers=True)
            prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pminus,yc,ts0,ts1,three_layers=True)
            prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pminus,yc,ts0,ts1,three_layers=True)
        else:
            prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pplus,yc)
            prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pplus,yc)
            prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pminus,yc)
            prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pminus,yc)
    
        prism_m_plus = prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,pplus,yc)
        prism_m_minus = prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,pminus,yc)
        gplus = g_function(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),gzw,prism_s_plus,prism_c_plus,prism_m_plus)
        gminus = g_function(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),gzw,prism_s_minus,prism_c_minus,prism_m_minus)

        pplus = np.copy(p)
        pminus = np.copy(p)
        
        for i in range(n):
            
            G[i,k] = (gplus[i] - gminus[i])/(2*dp)

    return G

def G_matrix_function_all(xmax,xmin,dy,edge,dp1,dp2,S0,dw,ds,dm,dref,dc,tw,p,yc,area,ipar,ts0=None,ts1=None,two_layers=False,three_layers=False,display=True):
    '''
    This function create the sensibility matrix considering the contribution 
    of all the prisms at a point of observation.
        
    Input
    xmax: float - maximum coordinate of the columns forming the
    interpretation model (North direction).
    xmin: float - minimum coordinate of the columns forming the
    interpretation model (North direction).
    dy: float - Thickness of columns forming the
    interpretation model (East direction).
    edge: float - Edge extension.
    dp1: float - ts and tm parameter variation to Finite Difference
    dp2: float - S0 and dS0 parameter variation to Finite Difference
    S0: numpy array 1D - isostatic compensation surface.
    dw: float - water density.
    ds: float - sediments density.
    dm: float - mantle density.
    dref: numpy array 1D - reference density.
    dc: numpy array 1D - crust density.
    tw: numpy array 1D - thickness of the water layer along
    the profile.
    p: numpy array 1D - vector parameters of the model.
    yc: numpy array 1D - coordinates of the center of the columns forming the
    interpretation model (East direction).
    area: numpy array 2D - model  limits = [ymin, ymax, zmax, zmin]
    ipar: integer - parameter indice for display of the finite difference step.
    ts0: numpy array 1D - if not None, thickness of the first sediment layer along
    the profile.
    ts1: numpy array 1D - if not None, thickness of the second sediment layer along
    the profile.
    two_layers: boolean - if true, the sediment layer has 2 layers. Otherwise it will has 1.
    three_layers: boolean - if true, the sediment layer has 3 layers. Otherwise it will has 1.
    display: boolean - if True, the function display the models of prisms with
    a parameter and its finite difference step.
        
    Output
    G: numpy array 2D - sensibility matriz.
    '''
    
    #initialization of variables
    x = np.zeros_like(yc)
    z = np.zeros_like(yc)-150.0
    n = len(yc)
    dS0 = p[n+n]
    
    gplus = np.zeros_like(yc)
    gminus = np.zeros_like(yc)
    pplus = np.copy(p)
    pminus = np.copy(p)
    
    G = np.zeros((n,2*n+1))
    
    #error messages
    assert np.min(tw) > 0.0, \
        'Bathymetry must be placed below the zero coordinate'
    assert ipar <= (2*n), \
        'parameter indice must be less than length of the parameter vector'
    
    #function implementation
    prism_w = prism_w_function(xmax,xmin,dy,edge,dw,dref,tw,yc)
    gzw = prism.gz(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),prism_w)
    
    if two_layers:
        tc = S0 - tw - ts0 - p[0:n] - p[n:n+n]
    elif three_layers:
        tc = S0 - tw - ts0 - ts1 - p[0:n] - p[n:n+n]
    else:
        tc = S0 - tw - p[0:n] - p[n:n+n]
    
    for k in range(2*n):
        dp = dp1
        if tc[k%n] <= dp:
            if k < n:
                pplus[k] = p[k]
                pminus[k] = p[k]-dp
            else:
                pplus[k] = p[k]+dp
                pminus[k] = p[k]   
        
        elif p[k] <= dp:
            if k < n:
                pplus[k] = p[k]+dp
                pminus[k] = p[k]
            else:
                pplus[k] = p[k]
                pminus[k] = p[k]-dp   
       
        else:
            pplus[k] = p[k]+dp
            pminus[k] = p[k]-dp 
               
        if two_layers:
            prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pplus,yc,ts0,two_layers=True)
            prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pplus,yc,ts0,two_layers=True)
            prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pminus,yc,ts0,two_layers=True)
            prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pminus,yc,ts0,two_layers=True)
        elif three_layers:
            prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pplus,yc,ts0,ts1,three_layers=True)
            prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pplus,yc,ts0,ts1,three_layers=True)
            prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pminus,yc,ts0,ts1,three_layers=True)
            prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pminus,yc,ts0,ts1,three_layers=True)
        else:
            prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pplus,yc)
            prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pplus,yc)
            prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pminus,yc)
            prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pminus,yc)

        prism_m_plus = prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,pplus,yc)
        prism_m_minus = prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,pminus,yc)
        
        gplus = g_function(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),gzw,prism_s_plus,prism_c_plus,prism_m_plus)
        gminus = g_function(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),gzw,prism_s_minus,prism_c_minus,prism_m_minus)

        pplus = np.copy(p)
        pminus = np.copy(p)
    
        for i in range(n):
    
            G[i,k] = (gplus[i] - gminus[i])/(2*dp)
    
    # Para os parametros S0 e dS0:
    for k in range(2*n,2*n+1):
        dp = dp2
        pplus[k] = p[k] + dp
        pminus[k] = p[k] - dp
                  
        if two_layers:
            prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pplus,yc,ts0,two_layers=True)
            prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pplus,yc,ts0,two_layers=True)
            prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pminus,yc,ts0,two_layers=True)
            prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pminus,yc,ts0,two_layers=True)
        elif three_layers:
            prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pplus,yc,ts0,ts1,three_layers=True)
            prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pplus,yc,ts0,ts1,three_layers=True)
            prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pminus,yc,ts0,ts1,three_layers=True)
            prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pminus,yc,ts0,ts1,three_layers=True)
        else:
            prism_s_plus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pplus,yc)
            prism_c_plus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pplus,yc)
            prism_s_minus = prism_s_function(xmax,xmin,dy,edge,ds,dref,tw,pminus,yc)
            prism_c_minus = prism_c_function(xmax,xmin,dy,edge,S0,dref,dc,tw,pminus,yc)
    
    
        prism_m_plus = prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,pplus,yc)
        prism_m_minus = prism_m_function(xmax,xmin,dy,edge,S0,dref,dm,pminus,yc)
        
        gplus = g_function(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),gzw,prism_s_plus,prism_c_plus,prism_m_plus)
        gminus = g_function(np.reshape(x,(n,)),np.reshape(yc,(n,)),np.reshape(z,(n,)),gzw,prism_s_minus,prism_c_minus,prism_m_minus)

        pplus = np.copy(p)
        pminus = np.copy(p)
        
        for i in range(n):
            G[i,k] = (gplus[i] - gminus[i])/(2*dp)  
    
    if display:
        yli = yc[ipar%n] - (0.5*dy)
        yri = yc[ipar%n] + (0.5*dy)
        if ipar < n:
            if two_layers:
                surfu = tw[ipar] + ts0[ipar] + (p[ipar] - dp)
                surfd = tw[ipar] + ts0[ipar] + (p[ipar] + dp)
            elif three_layers:
                surfu = tw[ipar] + ts0[ipar] + ts1[ipar] + (p[ipar] - dp)
                surfd = tw[ipar] + ts0[ipar] + ts1[ipar] + (p[ipar] + dp)
            else:
                surfu = tw[ipar] + (p[ipar] - dp)
                surfd = tw[ipar] + (p[ipar] + dp)
            dS0u = p[n+n]
            dS0d = p[n+n]
        elif ipar >= n and ipar < 2*n:
            surfu = S0 - (p[ipar] + dp)
            surfd = S0 - (p[ipar] - dp)
            dS0u = p[n+n]
            dS0d = p[n+n]
        else:
            yli = yc[0] - (0.5*dy)
            yri = yc[n-1] + (0.5*dy)
            #surfu = p[ipar-1] + p[ipar] - dp
            #surfd = p[ipar-1] + p[ipar] + dp
            surfu = S0 + p[ipar] - dp
            surfd = S0 + p[ipar] + dp
            
            dS0u = p[n+n]
            dS0d = p[n+n]
        if two_layers:
            polygons_display(surfu,surfd,dS0u,dS0d,yli,yri,dy,S0,dw,ds,dm,dref,dc,tw,p,yc,area,ts0,two_layers=True)
        elif three_layers:
            polygons_display(surfu,surfd,dS0u,dS0d,yli,yri,dy,S0,dw,ds,dm,dref,dc,tw,p,yc,area,ts0,ts1,three_layers=True)
        else:
            polygons_display(surfu,surfd,dS0u,dS0d,yli,yri,dy,S0,dw,ds,dm,dref,dc,tw,p,yc,area)

    return G
