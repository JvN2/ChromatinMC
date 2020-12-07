# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 10:18:02 2017

@author: noort
"""
import matplotlib as mpl

mpl.use(u'TkAgg')
mpl.interactive(True)

import numpy as np
import os
import subprocess
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from scipy import misc
# ChromatinMC modules:
import FileIO as fileio


POVRAY_BINARY = ("C:\Program Files\POV-Ray\\v3.7\\bin\\pvengine64.exe"
                 if os.name == 'nt' else "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64")


def init(rotation=[0, 0], plt_width=250, aspect_ratio=2):
    """
    Create Header for POVray file
    
    Parameters
    -------------
    rotation: ndarray of (2)
        [lattitude, longitude]
    lookat: ndarray of (3)
        focus of image
    camera: ndarray of (3)
        focus of image
    Returns
    ---------
    """

    HeaderTxt = """
// File: *.pov
// Creator: van Noort (c) 2017
// Version: POV-Ray Version 3

// Camera, Use aspectratio 1:#aspect_ratio#
camera {
    orthographic
    location <5000, 0, 0>
    look_at < 0, 0, 0>
    right <#right#>
    up <#up#>
}

// Light
light_source {< 15000, 0, 500> color rgb <1.0, 1.0, 1.0>}

background { color rgb <1, 1, 1> }

union {   
// Objects
   // union finish
   // change "T" to a higher value (up to 1.0, e.g. 0.8) and
   // uncomment "hollow" to make your object transparent
   #declare F = finish { specular 0.4 roughness 0.005
                         diffuse 0.8
                         ambient 0.2 }
   #declare T = 0;


// Insert objects here


   // no_shadow
   // hollow
   scale <1, -1, 1>
   rotate < -90, 0, 0 >
   translate <#translation#>

}
"""
    HeaderTxt = HeaderTxt.replace('#aspect_ratio#', '%d' % (aspect_ratio))
    HeaderTxt = HeaderTxt.replace('#rotation#', '%d, %d' % (rotation[0], rotation[1]))
    HeaderTxt = HeaderTxt.replace('#right#', '%d, %d, %d' % (plt_width, 0, 0))
    HeaderTxt = HeaderTxt.replace('#up#', '%d, %d, %d' % (0, plt_width * aspect_ratio, 0))
    HeaderTxt = HeaderTxt.replace('#translation#', '%d, %d, %d' % (0, 0, 0))
    return HeaderTxt


def add_sphere(pov_image, coord, radius=4, color=[1, 0, 0], transperancy=0.0):
    POVTxt = """
object {
   sphere { < #xyz# > #R# }
   texture {
      pigment {color rgbt <#rgb#, #t#>}
      finish {F}
           }
       }
"""

    rgb = np.array(color)

    POVTxt = POVTxt.replace('#xyz#', '%f, %f, %f' % (coord[0], coord[1], coord[2]))
    POVTxt = POVTxt.replace('#R#', '%f' % (radius))
    POVTxt = POVTxt.replace('#t#', '%f' % (transperancy))
    POVTxt = POVTxt.replace('#rgb#', '%f, %f, %f' % (rgb[0], rgb[1], rgb[2]))
    pov_image = pov_image.replace('// Insert objects here', '// Insert objects here %s' % POVTxt)
    return pov_image


def save(POVTxt, filename):
    with open(filename, 'w+') as f:
        f.write(POVTxt)
    return filename


def render(filename, height=None, width=None, quality=None, antialiasing=None, delete_pov=True):
    """ Renders a pov file with POV-Ray.
    
    Parameters
    ------------

    filename
      Name of the POV file for the output.
    
    height
      height in pixels

    width
      width in pixels

    """
    f = os.path.normpath(filename)
    cmd = [POVRAY_BINARY, filename]
    if height is not None: cmd.append('+H%d' % height)
    if width is not None: cmd.append('+W%d' % width)
    if quality is not None: cmd.append('+Q%d' % quality)
    if antialiasing is not None: cmd.append('+A%f' % antialiasing)
    cmd.append("Output_File_Type=N")
    cmd.append("/exit")
    cmd.append('-d')

    #    print cmd
    SW_HIDE = 0
    info = subprocess.STARTUPINFO()
    info.dwFlags = subprocess.STARTF_USESHOWWINDOW
    info.wShowWindow = SW_HIDE
    process = subprocess.Popen(cmd, stderr=subprocess.PIPE,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE, startupinfo=info)
    string = ''
    out, err = process.communicate(string.encode('ascii'))

    if process.returncode:
        raise IOError("POVRay rendering failed with the following error: " + err)

    if delete_pov:
        os.remove(filename)

    filename = fileio.change_extension(filename, 'png')
    return mpimg.imread(filename)


def show(filename):
    # plt.close()
    filename = fileio.change_extension(filename, 'png')
    img = misc.imread(filename)
    plt.imshow(img)
    plt.axis('off')
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    plt.draw()
    plt.show()
    # plt.pause(0.05)

    return

def main(filename, coord, colors, radius=10, range_A=[1000, 1000], offset_A=[0, 0, 0], width_pix=500,
                   showt=False):

    filename = fileio.change_extension(filename, 'pov')
    aspect_ratio = range_A[1] / float(range_A[0])
    pov_image = init(plt_width=range_A[0], aspect_ratio=aspect_ratio)
    # i = 0
    # j = 0
    offset = np.asarray(offset_A) - np.asarray((0, 0, range_A[1] / 2.0))
    for i, c in enumerate(coord):
        pov_image = add_sphere(pov_image, c + offset, color=colors[i], radius=radius)
    # for coord, t in zip(coord, transparency):
    #     if (i > len(colors) - 1):
    #         i = 0
    #     if (j > len(radius) - 1):
    #         j = 0
    #     for sphere in coord:
    #         pov_image = add_sphere(pov_image, sphere + offset, color=colors[i], radius=radius[j],
    #                                    transperancy=t)
    #     i += 1
    #     j += 1
    save(pov_image, filename=filename)
    render(filename, height=width_pix * aspect_ratio, width=width_pix)
    if showt:
        show(filename)
    filename = fileio.change_extension(filename, 'png')

    return filename

if __name__ == '__main__':

    # pov = init(plt_width=400, aspect_ratio=2)
    # r = 3
    # for i in range(100):
    #     pov = add_sphere(pov, [i, 0, 0], color='r', radius=r)
    #     pov = add_sphere(pov, [0, i, 0], color='g', radius=r)
    #     pov = add_sphere(pov, [0, 0, 2 * i], color='b', radius=r)
    filename = fileio.get_filename(root='test', ext='pov', incr=True)
    coord = [[10,0,0],[0,10,0],[0,0,10]]
    color = [[0.2,0,0],[0.4,0,0],[0.8,0,0]]
    print(main(filename, coord, colors=color, showt=True))
    # save(pov, filename)
    # size = 1000
    # render(filename, height=size, width=size)
    # show(filename)
