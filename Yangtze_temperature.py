# coding = 'utf-8'
from lytools import *
T = Tools()

this_root = 'C:\\Users\\leeya\\PycharmProjects\\ANUSPLIN\\Yangtze\\'
program_dir = join(this_root,'anusplin_program')
confdir = join(this_root,'conf\\template')

DEM_path = join(this_root,'conf\\chinadem0.0625.dem')
Yangtze_DEM_tif_path = join(this_root,'conf\\Yangtzedem0.0625.tif')
format_GCM_output_folder = join(this_root,'formated_GCM')
output_folder = join(this_root,'output')

def format_GCM():
    fpath = join(this_root, 'data\\GCM_data1.txt')
    T.mkdir(format_GCM_output_folder)
    extracted_DEM_dict = extract_DEM()
    f = open(fpath, 'r')
    lines = f.readlines()
    # f.close()
    dates = ''
    for line in lines:
        line = line.split('\n')[0]
        titles = line.split('\t')

        dates = titles[3:]
        break

    for i,date_str in tqdm(enumerate(dates),total=len(dates)):
        mon,year = date_str.split('/')
        date = year+mon
        outf = join(format_GCM_output_folder,date+'.dat')
        fw = open(outf,'w')
        for line in lines[1:]:
            vals_list = line.split('\t')
            # print(vals_list);exit()
            Station_ID = vals_list[0]
            lon = vals_list[1]
            lat = vals_list[2]
            val = vals_list[i+3]

            lon = float(lon)
            lat = float(lat)
            height = extracted_DEM_dict[Station_ID]
            val = float(val)
            # line_w = ' '.join([Station_ID,lon,lat,val])
            text = ' {:>5d} {:>7.4f} {:>6.4f}{:>8.2f}{:>8.2f}\n'.format(int(Station_ID), lon, lat, height, val)
            fw.write(text)

def gen_conf(date):

    for i in [1,2,3]:
        output_folder_date = join(output_folder,date)
        T.mkdir(output_folder_date,force=True)
        conf_path = join(confdir,f'{i}.conf')
        conf_path_r = open(conf_path,'r')
        lines = conf_path_r.readlines()
        conf_path_r.close()

        fw = open(join(output_folder_date,f'{i}.conf'),'w')
        for line in lines:
            line = line.split('\n')[0]
            new_line = line.replace('datadir',join(format_GCM_output_folder,date))
            new_line = new_line.replace('tempdir',join(output_folder_date,date))
            new_line = new_line.replace('dem_path',DEM_path)
            fw.write(new_line+'\n')
        fw.close()
        pass


def interpolate(date):
    output_folder_date = join(output_folder, date) + '\\'
    selnot = program_dir+'\\selnot.exe'
    splinb = program_dir+'\\splinb.exe'
    lapgrd = program_dir+'\\lapgrd.exe'
    conf1 = output_folder_date+'1.conf'
    conf2 = output_folder_date+'2.conf'
    conf3 = output_folder_date+'3.conf'
    cmd1 = selnot+' <{}> {}1.log'.format(conf1,output_folder_date+date)
    cmd2 = splinb+' <{}> {}2.log'.format(conf2,output_folder_date+date)
    cmd3 = lapgrd+' <{}> {}3.log'.format(conf3,output_folder_date+date)
    # print(cmd1)
    # print(cmd2)
    # print(cmd3)

    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    pass

def extract_DEM():
    D = DIC_and_TIF(tif_template=Yangtze_DEM_tif_path)
    spatial_dict = D.spatial_tif_to_dic(Yangtze_DEM_tif_path)

    fpath = join(this_root, 'data\\GCM_data1.txt')
    T.mkdir(format_GCM_output_folder)

    f = open(fpath, 'r')
    lines = f.readlines()
    f.close()
    sta_id_list = []
    lon_list = []
    lat_list = []
    for line in lines[1:]:
        vals_list = line.split('\t')
        # print(vals_list);exit()
        Station_ID = vals_list[0]
        sta_id_list.append(Station_ID)
        lon = vals_list[1]
        lat = vals_list[2]
        lon = float(lon)
        lat = float(lat)
        lon_list.append(lon)
        lat_list.append(lat)
    pix_list = D.lon_lat_to_pix(lon_list,lat_list)

    extracted_val_dict = {}
    for i,pix in enumerate(pix_list):
        val = spatial_dict[pix]
        sta = sta_id_list[i]
        extracted_val_dict[sta] = val

    return extracted_val_dict


def grd_to_tif(date):
    output_folder_date = join(output_folder, date)
    grd_fpath = join(output_folder_date,f'{date}.grd')
    outf = join(output_folder_date,f'{date}.tif')
    fr = open(grd_fpath).readlines()
    matrix = []
    for line in fr:
        line = line.split('\n')[0]
        vals = line.split()
        if len(vals) < 3:
            # print(vals)
            continue
        vals_float = [float(i) for i in vals]
        matrix.append(vals_float)
    matrix = np.array(matrix)
    # print(np.shape(matrix));
    # exit()
    longitude_start = 70
    latitude_start = 60
    pixelWidth = 0.0625
    pixelHeight = -0.0625

    Yangtze_DEM_tif_path_China = join(this_root, 'conf\\Yangtzedem0.0625_China.tif')
    Yangtze_DEM_arr = ToRaster().raster2array(Yangtze_DEM_tif_path_China)[0]
    mask_arr = Yangtze_DEM_arr < 0
    matrix[mask_arr] = np.nan
    array2raster(outf,longitude_start, latitude_start, pixelWidth, pixelHeight, matrix)

def array2raster(newRasterfn, longitude_start, latitude_start, pixelWidth, pixelHeight, array, ndv=-999999):
    cols = array.shape[1]
    rows = array.shape[0]
    originX = longitude_start
    originY = latitude_start
    # open geotiff
    driver = gdal.GetDriverByName('GTiff')
    if os.path.exists(newRasterfn):
        os.remove(newRasterfn)
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32)
    # outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_UInt16)
    # ndv = 255
    # Add Color Table
    # outRaster.GetRasterBand(1).SetRasterColorTable(ct)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    # Write Date to geotiff
    outband = outRaster.GetRasterBand(1)

    outband.SetNoDataValue(ndv)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    # outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(wkt_84())
    # Close Geotiff
    outband.FlushCache()
    del outRaster

def wkt_84():
    wkt_str = '''GEOGCRS["WGS 84",
ENSEMBLE["World Geodetic System 1984 ensemble",
    MEMBER["World Geodetic System 1984 (Transit)"],
    MEMBER["World Geodetic System 1984 (G730)"],
    MEMBER["World Geodetic System 1984 (G873)"],
    MEMBER["World Geodetic System 1984 (G1150)"],
    MEMBER["World Geodetic System 1984 (G1674)"],
    MEMBER["World Geodetic System 1984 (G1762)"],
    MEMBER["World Geodetic System 1984 (G2139)"],
    ELLIPSOID["WGS 84",6378137,298.257223563,
        LENGTHUNIT["metre",1]],
    ENSEMBLEACCURACY[2.0]],
PRIMEM["Greenwich",0,
    ANGLEUNIT["degree",0.0174532925199433]],
CS[ellipsoidal,2],
    AXIS["geodetic latitude (Lat)",north,
        ORDER[1],
        ANGLEUNIT["degree",0.0174532925199433]],
    AXIS["geodetic longitude (Lon)",east,
        ORDER[2],
        ANGLEUNIT["degree",0.0174532925199433]],
USAGE[
    SCOPE["Horizontal component of 3D system."],
    AREA["World."],
    BBOX[-90,-180,90,180]],
ID["EPSG",4326]]'''
    return wkt_str

def unify_raster(in_tif, out_tif, ndv=-999999):
    '''
    Unify raster to the extend of (70 140 10 60)
    '''
    insert_value = ndv
    array, originX, originY, pixelWidth, pixelHeight = ToRaster().raster2array(in_tif)
    # insert values to row
    top_line_num = abs((60. - originY) / pixelHeight)
    bottom_line_num = abs((-10. + originY + pixelHeight * len(array)) / pixelHeight)
    top_line_num = int(round(top_line_num, 0))
    bottom_line_num = int(round(bottom_line_num, 0))
    nan_array_insert = np.ones_like(array[0]) * insert_value
    top_array_insert = []
    for i in range(top_line_num):
        top_array_insert.append(nan_array_insert)
    bottom_array_insert = []
    for i in range(bottom_line_num):
        bottom_array_insert.append(nan_array_insert)
    bottom_array_insert = np.array(bottom_array_insert)
    if len(top_array_insert) != 0:
        arr_temp = np.insert(array, obj=0, values=top_array_insert, axis=0)
    else:
        arr_temp = array
    if len(bottom_array_insert) != 0:
        array_unify_top_bottom = np.vstack((arr_temp, bottom_array_insert))
    else:
        array_unify_top_bottom = arr_temp

    # insert values to column
    left_line_num = abs((70. - originX) / pixelWidth)
    right_line_num = abs((140. - (originX + pixelWidth * len(array[0]))) / pixelWidth)
    left_line_num = int(round(left_line_num, 0))
    right_line_num = int(round(right_line_num, 0))
    left_array_insert = []
    right_array_insert = []
    for i in range(left_line_num):
        left_array_insert.append(insert_value)
    for i in range(right_line_num):
        right_array_insert.append(insert_value)

    array_unify_left_right = []
    for i in array_unify_top_bottom:
        if len(left_array_insert) != 0:
            arr_temp = np.insert(i, obj=0, values=left_array_insert, axis=0)
        else:
            arr_temp = i
        if len(right_array_insert) != 0:
            array_temp1 = np.hstack((arr_temp, right_array_insert))
        else:
            array_temp1 = arr_temp
        array_unify_left_right.append(array_temp1)
    array_unify_left_right = np.array(array_unify_left_right)
    newRasterfn = out_tif
    ToRaster().array2raster(newRasterfn, 70, 60, pixelWidth, pixelHeight, array_unify_left_right, ndv=ndv)


def extend_Yangtze_DEM_to_China():
    outtif = Yangtze_DEM_tif_path.replace('.tif','_China.tif')
    unify_raster(Yangtze_DEM_tif_path,outtif,ndv=-32768)
    pass

def main():
    extend_Yangtze_DEM_to_China()
    format_GCM()

    datelist = []
    for f in os.listdir(format_GCM_output_folder):
        datelist.append(f.split('.')[0])
    for date in tqdm(datelist):
        gen_conf(date)
        interpolate(date)
        grd_to_tif(date)
        # exit()
    pass

if __name__ == '__main__':
    main()