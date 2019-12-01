# coding=gbk

import os
from tqdm import tqdm
import multiprocessing
from multiprocessing.pool import ThreadPool as TPool
import copy_reg
import types
this_root = r'D:/ANU_interpolation/'



class MUTIPROCESS:
    '''
    可对类内的函数进行多进程并行
    由于GIL，多线程无法跑满CPU，对于不占用CPU的计算函数可用多线程
    并行计算加入进度条
    '''
    def __init__(self,func,params):
        self.func = func
        self.params = params
        copy_reg.pickle(types.MethodType, self._pickle_method)
        pass

    def _pickle_method(self,m):
        if m.im_self is None:
            return getattr, (m.im_class, m.im_func.func_name)
        else:
            return getattr, (m.im_self, m.im_func.func_name)


    def run(self,process=6,process_or_thread='p',**kwargs):
        '''
        # 并行计算加进度条
        :param func: input a kenel_function
        :param params: para1,para2,para3... = params
        :param process: number of cpu
        :param thread_or_process: multi-thread or multi-process,'p' or 't'
        :param kwargs: tqdm kwargs
        :return:
        '''
        if 'text' in kwargs:
            kwargs['desc'] = kwargs['text']
            del kwargs['text']

        if process_or_thread == 'p':
            pool = multiprocessing.Pool(process)
        elif process_or_thread == 't':
            pool = TPool(process)
        else:
            raise IOError('process_or_thread key error, input keyword such as "p" or "t"')

        results = list(tqdm(pool.imap(self.func, self.params), total=len(self.params),**kwargs))
        pool.close()
        pool.join()
        return results






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
    # subprocess.call(cmd1)
    # subprocess.call(cmd2)
    # subprocess.call(cmd3)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    pass

def move_result_delete_temp_folder(mode,date):



    pass



def delete_0():
    # 删除降水小于0的像素

    pass

def kernel_run(params):
    mode, d = params
    gen_conf(mode, d)
    interpolate(mode, d)


def run():
    mode = 'pre'
    datelist = []
    fdir = this_root+'/data/daily/pre/'
    for f in os.listdir(fdir):
        datelist.append(f.split('.')[0])
    params = []
    for d in datelist:
        params.append([mode, d])
    MUTIPROCESS(kernel_run,params).run()


def main():

    run()
    pass


if __name__ == '__main__':
    main()