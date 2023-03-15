import pandas as pd
from tqdm import tqdm

def is_slf_intrsct(markup):    
    '''
    Функция is_slf_intrsct() определяет факт наличия пересечений
    сегментов внутри разметки markup. Требует один аргумент - разметку .deb.
    Возвращает True или False.
    '''
    
    # later should adapt for using lesser memory
    markup = markup.sort_values(by=[1], ascending=True)
    begin_list = list(markup[1])
    end_list = list(markup[2])
    slf_intrsct = False
    for zone_ in tqdm(range(1,len(begin_list))):
        if slf_intrsct == True:
            break
        if end_list[zone_ - 1] > begin_list[zone_]:
            slf_intrsct = True
    return slf_intrsct

def are_zones_intrsct(segm1, segm2):   
    '''
    Функция are_zones_intrsct() проверяет, пересекаются ли указанные 
    в аргументах отрезки zone1 и zone2. Отрезки должны списками или кортежами длины 2.
    Возвращает True или False.
    '''
    
    begin = 0
    end = 1
    intrsct = False
    
    if max(segm1[begin], segm2[begin]) < min(segm1[end], segm2[end]):
        intrsct = True
    
    return intrsct

def get_intrscts(markup1, markup2, ask_if=False):    
    '''
    Функция get_intrscts() строит матричное представление графа, построенного на
    отношениях пересечений между отрезками разметок markup1 и markup2.
    Возвращает объект pd.DataFrame. Если аргумент ask_if=True, то возвращает только True или False.
    '''
    
    begin_list_1 = list(markup1[1])
    end_list_1 = list(markup1[2])
    begin_list_2 = list(markup2[1])
    end_list_2 = list(markup2[2])
    df_intrsct = pd.DataFrame(0, index=[i for i in range(1,len(begin_list_2)+1)],
                               columns=[i for i in range(1,len(begin_list_1)+1)])
    for zone1 in tqdm(range(len(begin_list_1))):
        for zone2 in range(len(begin_list_2)):
            if are_zones_intrsct([begin_list_1[zone1], end_list_1[zone1]],
                                 [begin_list_2[zone2], end_list_2[zone2]]):
                if ask_if:
                    return True
                df_intrsct.loc[zone2+1, zone1+1] = 1
    return df_intrsct

def get_comm_intrsct(segm1, segm2):    
    '''
    Функция get_comm_intrsct() возвращает отрезок, являющийся пересечением
    отрезков segm1 и segm2.
    '''
    
    begin = 0
    end = 1
    
    if not are_zones_intrsct(segm1, segm2):
        print('Segments do not intersect')
        return 0
    
    intrsct_segm = [max(segm1[begin], segm2[begin]),
                    min(segm1[end], segm2[end])]
    
    return intrsct_segm
    
def build_intrscts(markup1, markup2, progr_off=False):   
    '''
    Функция build_intrscts() строит .deb разметку на основе пересечения двух разметок
    markup1 и markup2. Возвращает объект pd.DataFrame.
    '''
    
    begin_list_1 = list(markup1[1])
    end_list_1 = list(markup1[2])
    begin_list_2 = list(markup2[1])
    end_list_2 = list(markup2[2])
    markup1_len = 0
    markup2_len = 0
    intrsct_len = 0
    
    intrsct_markup = pd.DataFrame(columns=[0, 1, 2])
    
    # if markup2[0][0] != markup1[0][0]:
    if markup2.iat[0, 0] != markup1.iat[0, 0]:
        if input('Different chromosomes if files. If no,\
        process will be stopped, else will be used name from 1-st file') == 'no':
            return 0
    # chr_name = markup1[0][0]
    chr_name = markup1.iat[0, 0]
    
    for zone1 in tqdm(range(len(begin_list_1)), disable=progr_off):
        for zone2 in range(len(begin_list_2)):
            segm1 = [begin_list_1[zone1], end_list_1[zone1]]
            segm2 = [begin_list_2[zone2], end_list_2[zone2]]
            markup1_len += segm1[1] - segm1[0]
            markup2_len += segm2[1] - segm2[0]
            if are_zones_intrsct(segm1, segm2):
                intrsct_zone = get_comm_intrsct(segm1, segm2)
                intrsct_markup = pd.concat([intrsct_markup,
                                            pd.DataFrame([[chr_name] + intrsct_zone],
                                                         columns=[0, 1, 2])], ignore_index=True)
                intrsct_len += intrsct_zone[1] - intrsct_zone[0]
    intrsct_markup.sort_values(by=[1], ascending=True, inplace=True)
#     print('Общая разметка составляет ', 100*intrsct_len/markup1_len,'% и ',
#           100*intrsct_len/markup2_len, '% от исходных')
    return intrsct_markup

def segmentation_intrsct(markup1, markup2, chr_length=0, n_segments=10, progr_off=False):
    '''
    Функция segmentation_intrsct() строит .deb разметку на основе пересечения двух разметок
    markup1 и markup2, задействуя алгоритм деления на отрезки. Опционально использует
    длину хромосомы и количество отрезков.
    '''
    
    if not chr_length:
        chr_length = max(list(markup1[2])[-1], list(markup2[2])[-1])
    
    if (n_segments <= 0) or (n_segments != int(n_segments)):
        print('Неверное количество сегментов n, должно быть целое положительное')
        return 0
    
    intrsct_markup = pd.DataFrame(columns=[0, 1, 2])
    segments_lst = []
    for _s in range(1,n_segments+1):
        segments_lst.append(_s * (chr_length // n_segments))
    segments_lst[-1] += chr_length % n_segments
    for _s in tqdm(segments_lst, disable=progr_off):
        markup1_tmp = markup1.loc[(markup1[2] <= _s)]
        markup2_tmp = markup2.loc[(markup2[2] <= _s)]
        markup1_border = markup1.loc[(markup1[2] > _s) & (markup1[1] < _s)]
        markup2_border = markup2.loc[(markup2[2] > _s) & (markup2[1] < _s)]
        if markup1_tmp.empty or markup2_tmp.empty:
            continue
        if not markup1_border.empty:
            intrsct_markup = pd.concat([intrsct_markup, build_intrscts(markup1_border, markup2_tmp, progr_off=True)])
        if not markup2_border.empty:
            intrsct_markup = pd.concat([intrsct_markup, build_intrscts(markup1_tmp, markup2_border, progr_off=True)])
        markup1 = markup1.loc[(markup1[2] > _s)]
        markup2 = markup2.loc[(markup2[2] > _s)]
        intrsct_markup = pd.concat([intrsct_markup, build_intrscts(markup1_tmp, markup2_tmp, progr_off=True)])
    
    return intrsct_markup