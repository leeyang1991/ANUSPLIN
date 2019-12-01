# coding=gbk

import os
from tqdm import tqdm

this_root = r'D:/ANU_interpolation/'


def format_pre():
    fdir = this_root+'data/meteo/PRE18/'
    out_dir = this_root+'data/daily/pre/'
    for f in tqdm(os.listdir(fdir)):
        fr = open(fdir+f,'r')
        lines = fr.readlines()
        split_daily_pre(lines,out_dir)
    pass


def split_daily_pre(lines,out_dir):
    date_dic = {}
    for line in lines:
        line = line.split('\n')[0]
        line_split = line.split()
        sta, lat, lon, height, year, mon, day, _, _, val, _, _, _ = line_split
        date = '{}{}{}'.format(year, '%02d' % int(mon), '%02d' % int(day))
        date_dic[date] = []

    for line in lines:
        line = line.split('\n')[0]
        line_split = line.split()
        sta, lat, lon, height, year, mon, day, _, _, val, _, _, _ = line_split
        if float(val) > 30000:
            continue
        lat = float(lat)/100.
        lon = float(lon)/100.
        height = float(height)/10.
        val = float(val)/10.
        date = '{}{}{}'.format(year,'%02d'%int(mon),'%02d'%int(day))
        date_dic[date].append([sta, lon, lat, height, val])
        # print text
    for date in date_dic:
        f = open(out_dir+date+'.dat','w')
        for i in date_dic[date]:
            sta, lon, lat, height, val = i
            text = ' {} {}\t{}\t{}\t{}\n'.format(sta, '%0.2f'%lon, '%0.2f'%lon, '%0.1f'%height, val)
            f.write(text)
        f.close()
    pass




def format_tmp():
    fdir = this_root+'data/meteo/TEM18/'
    out_dir = this_root+'data/daily/tmp/'
    for f in tqdm(os.listdir(fdir)):
        fr = open(fdir+f,'r')
        lines = fr.readlines()
        split_daily_pre(lines,out_dir)
    pass


def split_daily_tmp(lines,out_dir):
    date_dic = {}
    for line in lines:
        line = line.split('\n')[0]
        line_split = line.split()
        sta, lat, lon, height, year, mon, day, mean, max_v, min_v, _, _, _ = line_split
        date = '{}{}{}'.format(year, '%02d' % int(mon), '%02d' % int(day))
        date_dic[date] = []

    for line in lines:
        line = line.split('\n')[0]
        line_split = line.split()
        sta, lat, lon, height, year, mon, day, mean, max_v, min_v, _, _, _ = line_split
        if float(mean) > 30000:
            continue
        lat = float(lat)/100.
        lon = float(lon)/100.
        height = float(height)/10.
        mean = float(mean)/10.
        date = '{}{}{}'.format(year,'%02d'%int(mon),'%02d'%int(day))
        date_dic[date].append([sta, lon, lat, height, mean])
        # print text
    for date in date_dic:
        f = open(out_dir+date+'.dat','w')
        for i in date_dic[date]:
            sta, lon, lat, height, val = i
            text = ' {} {}\t{}\t{}\t{}\n'.format(sta, '%0.2f'%lon, '%0.2f'%lon, '%0.1f'%height, val)
            f.write(text)
        f.close()
    pass






def main():
    format_tmp()
    pass


if __name__ == '__main__':
    main()