# coding=gbk

import os
import shutil
import subprocess

this_root = r'D:/ANU_interpolation/'

def mkdir(fdir):
    if not os.path.isdir(fdir):
        os.makedirs(fdir)

def gen_conf(mode,date):
    data_dir = this_root+'data/daily/{}/'.format(mode)
    confdir = this_root+'conf/'
    for i in [1,2,3]:
        temp_folder = this_root+'temp_folder/'
        temp_folder_date = temp_folder+mode+'/'+date+'/'
        mkdir(temp_folder_date)
        # 1
        f1 = this_root+'conf/template/{}.conf'.format(i)
        fr1 = open(f1,'r')
        lines = fr1.readlines()
        fr1.close()
        fnew1 = open(temp_folder_date+'{}.conf'.format(i),'w')
        for line in lines:
            line = line.split('\n')[0]
            new_line = line.replace('datadir',data_dir+date)
            new_line = new_line.replace('tempdir',temp_folder_date+date)
            new_line = new_line.replace('confdir',confdir)
            fnew1.write(new_line+'\n')
        fnew1.close()
        pass

def interpolate(mode,date):
    program_dir = this_root+'anusplin_program/'
    temp_folder = this_root + 'temp_folder/'
    temp_folder_date = temp_folder + mode + '/' + date + '/'
    selnot = program_dir+'selnot.exe'
    splinb = program_dir+'splinb.exe'
    lapgrd = program_dir+'lapgrd.exe'
    conf1 = temp_folder_date+'1.conf'
    conf2 = temp_folder_date+'2.conf'
    conf3 = temp_folder_date+'3.conf'
    cmd1 = selnot+' <{}> {}.log'.format(conf1,temp_folder_date+date)
    cmd2 = splinb+' <{}> {}.log'.format(conf2,temp_folder_date+date)
    cmd3 = lapgrd+' <{}> {}.log'.format(conf3,temp_folder_date+date)

    pass

def move_result_delete_temp_folder(mode,date):



    pass



def delete_0():
    # 删除降水小于0的像素

    pass



def main():
    gen_conf('pre','19910618')
    interpolate('pre','19910618')
    pass


if __name__ == '__main__':
    main()