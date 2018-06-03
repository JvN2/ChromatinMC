import matplotlib as mpl

try:
    mpl.use(u'TkAgg')
    mpl.interactive(False)
except:
    pass

import time
import re
import os as os
import numpy as np
import pandas as pd
import sys
from collections import OrderedDict
from etaprogress.progress import ProgressBar
from lmfit import Parameters
from datetime import datetime
import getpass
import distutils.dir_util
import glob
import cv2
from matplotlib import pyplot as plt
from helixmc.pose import HelixPose
# ChromatinMC modules:
import NucleosomeMC as nMC
import POVutils as pov

default_folder = 'D:\\users\\'
kT = 41


def _merge(df1, df2):
    r = list(df2.index.values)[0]
    if r in list(df1.index.values):  # existing row
        df3 = df1
        for c in list(df2):
            df3.at[r, c] = df2.at[r, c]
    else:  # new row
        df3 = df1.append(df2)
    return df3


def report_progress(value, title='', init=False):
    global bar, start_time
    if init:
        start_time = time.time()
        print(datetime.now().strftime('>>> Start [{0}] @ %Y-%m-%d %H:%M:%S'.format(title)))
        bar = ProgressBar(value, max_width=80)
    else:
        bar.numerator = value
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))

    sys.stdout.write('\r>>> {0} = {1} {2}'.format(elapsed_time, bar, title))
    sys.stdout.flush()
    if bar.done:
        print ('\n')
    return


def get_filename(root=None, ext='dat', incr=False, sub=False, folder=False, list='', wildcard='*', date='today'):
    global working_directory, root_name, file_number, sub_number

    if 'root_name' not in globals():
        if root is None:
            root_name = 'data'
    if root is not None:
        root_name = root

    today = datetime.now().strftime('%Y%m%d')
    if date is 'today':
        folder_name = today
    elif date is 'yesterday':
        folder_name = (pd.to_datetime('Today') - pd.Timedelta('1 days')).strftime('%Y%m%d')
    else:
        foldername = date
    if 'working_directory' not in globals():
        user = getpass.getuser()
        working_directory = default_folder + '{0:s}\\data\\{1:s}\\'.format(user, today)
    filename = working_directory + wildcard + '.' + ext
    if date is not 'today':
        filename = default_folder + '{0:s}\\data\\{1:s}\\'.format(user, foldername) + wildcard + '.' + ext
    if list == 'all':
        list_of_files = glob.glob(filename)
        return list_of_files
    if list == 'last':
        list_of_files = glob.glob(filename)
        filename = max(list_of_files, key=os.path.getctime)
        return [filename]

    if 'file_number' not in globals():
        try:
            filename = working_directory + root_name + '_' + '*.*'
            file_number = int(re.split(r'_(\d\d\d)', (sorted(glob.glob(filename)))[-1])[-2])
        except:
            file_number = 0
    if 'sub_number' not in globals():
        sub_number = 0
    if incr:
        if sub:
            sub_number += 1
        else:
            file_number += 1
            sub_number = 0
    if sub:
        filename = working_directory \
                   + root_name + '_' + str(file_number).zfill(3) + '\\' \
                   + root_name + '_' + str(file_number).zfill(3) + '_' + str(sub_number).zfill(4) + '.' + ext
    else:
        filename = working_directory \
                   + root_name + '_' + str(file_number).zfill(3) + '.' + ext

    distutils.dir_util.mkpath(os.path.dirname(os.path.abspath(filename)))

    if folder:
        filename = os.path.dirname(filename)
    return filename


def change_extension(filename, extension):
    base, _ = os.path.splitext(filename)
    if extension.count('.') == 0:
        extension = '.' + extension
    return base + extension


def read_xlsx_row(filename, dataset, pars=None):
    # Read a row in Excel sheet
    sets, files, _ = contents_xlsx(filename)
    df = pd.read_excel(filename, 'Value')
    for index in df.index.values:
        if int(index.split(' > ')[0]) == dataset:
            row = index

    if pars is None:
        pars = Parameters()
        for par in list(df):
            pars.add(str(par), value=df.at[row, par])
        return pars, files[sets == dataset]

    for p in pars:
        try:
            pars[p].value = df.at[row, p]
        except Exception as e:
            pass
        df = pd.read_excel(filename, 'Min')
        try:
            pars[p].min = df.at[row, p]
        except:
            pass
        df = pd.read_excel(filename, 'Max')
        try:
            pars[p].max = df.at[row, p]
        except:
            pass
        df = pd.read_excel(filename, 'StErr')
        try:
            if df.at[row, p] == 0:
                pars[p].vary = False
            else:
                pars[p].vary = True
        except:
            pass
    return pars, files[sets == dataset]


def read_xlsx_collumn(filename, param):
    filename = change_extension(filename, 'xlsx')
    df = pd.read_excel(filename, 'Value')
    try:
        data = np.asarray(df[param])
        return data
    except:
        print ('Parameter [{:}] not found'.format(param))
    return


def write_xlsx_column(filename, param, data, report_file=None):
    if not (report_file):
        report_file = filename
    filename = change_extension(filename, 'xlsx')
    df = pd.read_excel(filename, 'Value')
    df[param] = data
    writer = pd.ExcelWriter(report_file.format('openpyxl'), engine='openpyxl')
    df.to_excel(writer, sheet_name=sheet)
    writer.close()
    return


def write_xlsx_row(filename, dataset, pars, report_file=None):
    dataset = int(np.clip([dataset],0,np.inf)[0])
    if not (report_file):
        report_file = filename
    report_file = change_extension(report_file, 'xlsx')
    writer = pd.ExcelWriter(report_file.format('openpyxl'), engine='openpyxl')
    df_index = str(dataset) + ' > ' + filename

    sheets = ['Value', 'StErr', 'Min', 'Max']
    for sheet in sheets:
        p = OrderedDict()
        if sheet == 'StErr':
            for k, v in pars.items():
                try:
                    p[k] = v.ster
                except:
                    if v.vary:
                        p[k] = np.inf
                    else:
                        p[k] = 0
        elif sheet == 'Min':
            for k, v in pars.items():
                p[k] = v.min
        elif sheet == 'Max':
            for k, v in pars.items():
                p[k] = v.max
        elif sheet == 'Value':
            for k, v in pars.items():
                p[k] = v.value

        df = pd.DataFrame(p, index=[df_index])
        try:
            old_dfs = pd.read_excel(report_file, sheet)
            df = _merge(old_dfs, df)
        except:
            pass
        df.to_excel(writer, sheet_name=sheet)

    writer.close()
    return df_index


def contents_xlsx(filename):
    filename = change_extension(filename, 'xlsx')
    df = pd.read_excel(filename, 'Value')
    datasets = []
    files = []
    for txt in list(df.index.values):
        datasets.append(int(txt.split(' > ')[0]))
        files.append(str(txt.split(' > ')[1]))
        par_names = []
        for item in list(df):
            par_names.append(str(item))
    return datasets, files, par_names


def create_movie(image_files, fps=5.0, delete=True, filename=None, reverse=True, titles=None):
    image_files = sorted(image_files)
    if reverse:
        image_files.append(reversed(image_files))
        titles.append(reversed(titles))

    if filename == None:
        filename = image_files[0]
    filename = change_extension(filename, 'avi')

    frame = cv2.imread(image_files[0])
    height, width, channels = frame.shape
    bottomLeftCornerOfText = (20, height - 20)
    font = cv2.FONT_HERSHEY_DUPLEX
    fontScale = 3
    fontColor = (0, 0, 0)
    lineType = 2

    fourcc = cv2.VideoWriter_fourcc(*'XVID')  # Be sure to use lower case
    out = cv2.VideoWriter(filename, fourcc, fps, (width, height))

    for i, image in enumerate(image_files):
        frame = cv2.imread(image)
        try:
            cv2.putText(frame, titles[i],
                        bottomLeftCornerOfText,
                        font,
                        fontScale,
                        fontColor,
                        lineType)
        except:
            pass

        out.write(frame)  # Write out frame to video

    if delete:
        for image in image_files:
            os.remove(image)
    out.release()
    print('>>> Movie saved as: {} \n'.format(filename))
    return filename


def save_plot(data, ax_labels=None, grid=None, xrange=None, yrange=None, save=True,
              filename=None, transpose=False):
    #  used format of data: [x, y, line, selected]
    if filename is None:
        filename = get_filename(incr=True)
    plt.close()
    plt.figure(figsize=(4, 3))
    plt.axes([0.15, 0.15, .8, .75])
    if not transpose:
        xsel = data[0][data[3] != 0]
        ysel = data[1][data[3] != 0]
        xunsel = data[0][data[3] == 0]
        yunsel = data[1][data[3] == 0]
        xline = data[0]
        yline = data[2]
    else:
        xsel = data[1][data[3] != 0]
        ysel = data[0][data[3] != 0]
        xunsel = data[1][data[3] == 0]
        yunsel = data[0][data[3] == 0]
        xline = data[2]
        yline = data[0]
    plt.scatter(xunsel, yunsel, s=10, facecolors='none', edgecolors='grey')
    plt.scatter(xsel, ysel, s=10, facecolors='none', edgecolors='b')
    plt.plot(xline, yline, color='k', linewidth=1.2)

    if grid is not None:
        for g in grid:
            if not transpose:
                plt.plot(xline, g, color='k', linestyle=':', linewidth=0.8)
            else:
                plt.plot(g, yline, color='k', linestyle=':', linewidth=0.8)

    if ax_labels is not None:
        plt.xlabel(ax_labels[0])
        plt.ylabel(ax_labels[1])
    if xrange is not None:
        plt.xlim(xrange)
    if yrange is not None:
        plt.ylim(yrange)
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.title(filename.split('\\')[-2] + '\\' + filename.split('\\')[-1].split('.')[0] + '\n', loc='left',
              fontdict={'fontsize': 8})
    # plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.draw()
    plt.pause(0.05)
    if save:
        filename = change_extension(filename, 'jpg')
        plt.savefig(filename, dpi=600, format='jpg')
    return filename


def create_pov(filename, coords, range_A=[1000, 1000], offset_A=[0, 0, 0], width_pix=500, colors=None, radius=None,
               show=False, transparency=None):
    if radius is None:
        radius = np.append([10], (np.ones(8) * 4))
        radius = np.append(radius, 3)
        radius = np.append(radius, 2)
    if colors is None:
        colors = 'kbbggryrycy'
    if transparency is None:
        transparency = np.zeros(len(coords))
    filename = change_extension(filename, 'pov')
    aspect_ratio = range_A[1] / float(range_A[0])
    pov_image = pov.init(plt_width=range_A[0], aspect_ratio=aspect_ratio)
    i = 0
    j = 0
    offset = np.asarray(offset_A) - np.asarray((0, 0, range_A[1] / 2.0))
    for coord, t in zip(coords, transparency):
        if (i > len(colors) - 1):
            i = 0
        if (j > len(radius) - 1):
            j = 0
        for sphere in coord:
            pov_image = pov.add_sphere(pov_image, sphere + offset, color=colors[i], radius=radius[j], transperancy=t)
        i += 1
        j += 1
    pov.save(pov_image, filename=filename)
    pov.render(filename, height=width_pix * aspect_ratio, width=width_pix)
    if show:
        pov.show(filename)
    filename = change_extension(filename, 'png')
    # print('File created: {:s}'.format(filename))
    return filename


def plot_dna(dna_pose1, origin_index=0, color='blue', update=False, title='', range_nm=100, save=False, wait=0,
             dna_pose2=None):
    global ax, fig, scale

    tf = nMC.get_transformation(nMC.get_of(dna_pose1, origin_index))
    coords = nMC.apply_transformation(dna_pose1.coords, tf) / 10.0

    if update:
        ax.clear()
        plt.title(title, loc='left')
    else:
        plt.close()
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(111, projection='3d')
        scale = range_nm / 2
        plt.title(title, loc='left')
        plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.3)
    pointsize = 1 * np.sqrt(scale)

    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    ax.set_zlabel('z (nm)')

    half = len(coords) / 2
    ax.scatter(coords[:, 0][:half], coords[:, 1][:half], coords[:, 2][:half], s=pointsize, c=color, alpha=0.5)
    ax.scatter(coords[:, 0][half:], coords[:, 1][half:], coords[:, 2][half:], s=pointsize, c='red', alpha=0.5)

    ax.scatter(coords[0, 0], coords[0, 1], coords[0, 2], s=pointsize, c='red', alpha=0.55)
    ax.scatter([0], [0], [0], s=pointsize * 5, c='k', alpha=0.55)

    if dna_pose2 is not None:
        coords = nMC.apply_transformation(dna_pose2.coords, tf) / 10.0
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=pointsize, c='black', alpha=0.25)

    plt.xlim((-scale, scale))
    plt.ylim((-scale, scale))
    ax.set_zlim(0, 2 * scale)

    plt.draw()
    plt.pause(0.05)
    if save:
        filename = get_filename(incr=True, sub=True, ext='png')
        plt.savefig(filename, dpi=600, format='png')
        print('plot saved in file: ', filename)
    time.sleep(wait)
    return


def create_pov_movie(filename, origin_frame=0, fps=5, reverse=False, frame=[], octamers=False, overwrite=False):
    pix = 1500
    radius = [10, 7, 7, 35]
    sets, filenames, _ = contents_xlsx(filename)

    for i, f in enumerate(filenames):
        filenames[i] = change_extension(f, 'png')

    if overwrite:
        existing_files = []
    else:
        path = (filename.split('.')[0])
        existing_files = glob.glob(path + '\\*.png')
    image_files = []

    names = ['n_nuc', 'NRL', 'dyad0']
    nucs = []
    for name in names:
        nucs.append(read_xlsx_collumn(filename, name))
    nuc = (np.asarray(nucs).T[0])

    nucl = nMC.NucPose()
    nucl.from_file('1KX5')
    n_of_coords_even = []
    octa_coords = []

    j = 0
    print('>>> {}'.format(filename))
    report_progress((len(filenames) - len(existing_files) - 1), title='create_pov_movie', init=True)
    for file in filenames:
        if file in existing_files:
            image_files.append(file)
        else:
            report_progress(j, title=file)
            try:
                dna = HelixPose.from_file(change_extension(file, 'npz'))
                origin_of = nMC.get_of(dna, origin_frame)
                tf = nMC.get_transformation(origin_of)
                coords = nMC.apply_transformation(dna.coords, tf)

                n_of_coords_odd = np.zeros((1, 3))
                n_of_coords_even = np.zeros((1, 3))
                if len(frame) is not 0:
                    for n in range(nuc[0]):
                        dyad = nuc[2] + n * nuc[1]
                        n_of = nMC.get_nuc_of(dna.coords, dna.frames, dyad, nucl)
                        if n % 2 is 0:
                            n_of_coords_even = np.concatenate((n_of_coords_even, nMC.of2axis(n_of, length=frame)))
                        else:
                            n_of_coords_odd = np.concatenate((n_of_coords_odd, nMC.of2axis(n_of, length=frame)))
                    n_of_coords_odd = nMC.apply_transformation(n_of_coords_odd[1:], tf)
                    n_of_coords_even = nMC.apply_transformation(n_of_coords_even[1:], tf)

                if octamers:
                    for n in range(nuc[0]):
                        dyad = nuc[2] + n * nuc[1]
                        n_of = nMC.get_nuc_of(dna.coords, dna.frames, dyad, nucl)
                        if n is 0:
                            octa_coords = [n_of[0]]
                        else:
                            octa_coords = np.concatenate((octa_coords, [n_of[0]]))
                    octa_coords = nMC.apply_transformation(octa_coords, tf)

                file = create_pov(file, [coords, n_of_coords_even, n_of_coords_odd, octa_coords], range_A=[1000, 2500],
                                  offset_A=[0, 0, 150],
                                  show=False, width_pix=pix, colors='kcyr', radius=radius)

                image_files.append(file)
                j += 1
            except:
                pass

    forces = read_xlsx_collumn(filename, 'F_pN')
    titles = []
    for force in forces:
        titles.append('F = {:2.1f} pN'.format(force))

    create_movie(image_files, fps=fps, delete=False, filename=filename, reverse=reverse,
                 titles=titles)
    return


def main():
    filename = 'E:\\users\\noort\\data\\20180321\\12nucs_003.dat'
    analyze_step_params(filename)
    return


if __name__ == "__main__":
    # execute only if run as a script
    main()
    # create_pov_mp4('E:\\Users\\Noort\\data\\20180315\\2nucs_001', origin_frame=0, reverse=False)
    # filename = get_filename(incr=True, base='2nucs')
    # print(filename)
