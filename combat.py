from cctbx import crystal
from cctbx import miller
from iotbx import scalepack
from cctbx.array_family import flex
import sys
import iotbx
from iotbx import pdb
import math
from cctbx.array_family import flex
from cctbx import xray
from iotbx.scalepack import merge
from libtbx.utils import Sorry
#from iotbx.scalepack import merge

def run(args):

  vars = get_input_dict(args)
  
  if vars['xia2'] == 'yes':
    unit_cell_params = xia2(vars['bragg'])

  else:
    unit_cell_params = [vars['cella'], vars['cellb'], vars['cellc'], vars['anglea'], vars['angleb'],vars['anglec'],vars['spacegroup']]

  files = diffuse_conversion(vars['diffuse'])


  #print '***************************************'
  #print files
  #print '***************************************'
  
  if vars['lunus'] == 'yes':
    os.system("mkdir " + vars['diffuse'] +"/processed")

    for file in files:
      frame_processing(vars['diffuse'],file,vars['punchim'],vars['thrshim'],vars['polarim'],vars['modeim'])

    for file in files:
      frame_averaging(vars['diffuse'] + '/processed/proc.' + file)

    print 'The reference frame is %s' %vars['reference']

    reference_frame(vars['diffuse'] + '/processed/proc.', vars['reference'])

    for file in files:
      print 'The file currently being processed is %s' %file
      frame_scaling(vars['diffuse'] + '/processed/proc.', file)

    genlat(vars['diffuse'],files)

  if vars['genlat'] == 'yes':

    args = args_generator(vars['diffuse'], vars['index_1'], vars['index_2'], vars['index_3'], unit_cell_params[0], unit_cell_params[1], unit_cell_params[2], vars['resolution'], vars['lattice_name'], vars['processors'],unit_cell_params, vars['known_setting'], vars['file_format'])

    print args
    
    os.system('libtbx.python genlat_labelit.py ' + args)

    #if vars['file_format'] == hkl:
       #data = dic_construct()

    #"args" are the arguments normally passed to integrate.py on the command line (put into get_input_dict command)
    #diffuse_integration(args)

  if vars['symmetry'] == 'yes':

    map_symmetry_extension(vars['lattice_name'], unit_cell_params)
  
    friedel_hkl(vars['lattice_name'])


  #What is the file name that will be passed to aniso_convert()? It's the variable "lattice_name" with "_raw" appended (less any changes due to symmetry extensions)
  if vars['anisotropic'] == 'yes':
    aniso_convert(vars['lattice_name'], vars['cella'], vars['cellb'], vars['cellc'], vars['resolution'], vars['file_format'])
    print vars['file_format'] #Will use second half of string addition to adjust for proper file name

  if vars['statistics'] == 'yes':
    os.system('phenix.reflection_statistics isotropic.mtz > isotropic_statistics.txt')
    os.system('phenix.reflection_statistics anisotropic.mtz > anisotropic_statistics.txt')
  
def get_input_dict(args):
# return dictionary of keywords and their values

  fin=open(args, 'r')

  lines = fin.readlines()

  input_dict={}
  for line in lines:
  #print arg   
    spl=line.split("=")
    #print spl
    if len(spl)==2:
      input_dict[spl[0]]=spl[1][:-1]
      print input_dict[spl[0]]

  input_dict['resolution'] = float(input_dict['resolution'])
  input_dict['processors'] = int(input_dict['processors'])
  input_dict['cella'] = float(input_dict['cella'])
  input_dict['cellb'] = float(input_dict['cellb'])
  input_dict['cellc'] = float(input_dict['cellc'])
  input_dict['anglea'] = int(input_dict['anglea'])
  input_dict['angleb'] = int(input_dict['angleb'])
  input_dict['anglec'] = int(input_dict['anglec'])
  input_dict['known_setting']=int(input_dict['known_setting'])

  return input_dict

def cbf2img(file):

  from dxtbx.format.Registry import Registry
  import sys, os

  f = None
  if file.split(".")[-1].lower() == "cbf" and os.path.exists(file):
    if f is None:
      f = Registry.find(file)
    img = f(file)
    db = img.get_detectorbase()
    db.readHeader()
    db.read()
    db.show_header()
    destpath = file.rstrip(".cbf") + ".img"
    print "Writing %s as %s"%(file,destpath)

    db.debug_write(destpath)

def xia2(filepath):
#Worry about changing directories?
    #os.system("cd %s"%(filepath))
    os.system("xia2 -3dii "+filepath)
    
    fin = open("xia2-summary.dat", 'r')
    lines = fin.readlines()
    unit_cell = lines[len(lines)-2]
    space_group = lines[len(lines)-1]

    q = str(unit_cell)
    a = q.split()
    del a[0]

    b = str(space_group)

    d = b.translate(None, 'Spacegroup:')
    d = d[:-1]
    d = d[1:]
    a.append(d)

    return a


def diffuse_conversion(location):
  dirListing = os.listdir(location)
  files = []
  #Better to use dictionary for this? -> Create key for each image and have corresponding
  for file in dirListing:
      if file.endswith('.cbf'):
          print file
          new_file = location + '/' + file
          print new_file
          cbf2img(new_file)
          files.append(file.rstrip(".cbf") + ".img")
      elif file.endswith('.img'):
          files.append(file)

  return files
        

def frame_processing(filepath,file,punch,thr,pol,mode):
#frame_processing is equivalent to proc.all
    
  p = filepath + '/'
  
  #os.system("mkdir " + p +"processed")

  punch_var=process(punch)

  thr_var=process(thr)

  pol_var=process(pol)

  mode_var=process(mode)

  print file

  print "punchim " + p + file + " " + punch_var + "image0.img"
    
  #punchim removes pixels within a specified XY region
  os.system("punchim " + p + file + " " + punch_var + p + "image0.img")

  #thrshim removes pixels above and below a given threshold
  os.system("thrshim " + p + "image0.img" + " " + thr_var + " " + p + "image.img")

  #polarim corrects for beam polarization
  os.system("polarim "+ p + "image.img " + p + "image00.img " + pol_var)


  #normim corrects for solid-angle normalization and detector-face rotation in a diffraction image
  os.system("normim " + p +  "image00.img " + p +  "image1.img")

  #modeim removes the Bragg peaks from an image by mode filtering using a specified mask size
  os.system("modeim "+ p +  "image1.img " + p + "image2.img " + mode_var)

  #Move the processed images to a new folder
  os.system("cp " + p +  "image2.img " + p +  "processed/proc." + file)

  #Clean up dummy images
  os.system("rm " + p + "image.img; rm " + p +  "image1.img; rm " +p +"image2.img; rm " + p + "image00.img")
  os.system("echo "+ file + " was successfully processed!")
  return

def process(string):
  data = string.split(',')

  line = ''

  for item in data:
    line += item
    line += ' '

  return line


def frame_averaging(frame):
#frame_averaging is equivalent to proc.avgim

  file_prefix = frame.rstrip('.img')
  file_rf = file_prefix + ".rf"
  file_asc = file_prefix + ".asc"
  file_s_img = file_prefix + "s.img"
  file_s_sq = file_prefix + 's.sqr.rf'
  
  #avgrim calculates the average intensity vs. radius for an input image
  os.system("avgrim "+ frame +" " + file_rf)
  os.system("echo "+ frame + " avgrim works!")

  #binasc converts the average intensity vs. radius binary file to ASCII format (???)
  os.system("binasc 2 < " + file_rf + " > " + file_asc)
  os.system("echo "+ frame + " binasc works!")

  #subrfim subtracts I(r) (in this case the average radial intensity) from the diffraction image
  os.system("subrfim " + file_rf + " " + frame +" " +file_s_img)
  os.system("echo "+ frame + " subrfim works!")
  #avsqrim calculates the average *squared* pixel value vs. radius
  os.system("avsqrim " + file_s_img  + " " + file_s_sq)
  os.system("echo "+ frame + " avsqrim works!")

  #binasc converts the average intensity *squared* vs. radius binary file as: float input -> ascii output w/ leading index
  os.system("binasc 2 < " + file_s_sq + " > " + file_prefix + "s.sqr.asc")
  os.system("echo "+ frame + " binasc works!")

  #Cleaning up dummy images
  os.system("rm " + file_s_sq+"; rm "+ file_s_img)
  os.system("echo "+ frame + " file cleanup works!")

def reference_frame(prefix,file):
#reference_frame is equivalent to proc.makeref
# Can write in new variables for proc.makeref arguments
  data = file.split('.')
  image = prefix + data[0] +'.asc'
  os.system('echo IMAGE = ' +image)
  os.system('./proc.makeref '+ image + ' 250 1000' )
  os.system("echo "+ image + " reference frame creation works!")


def frame_scaling(prefix,file):

  data = file.split('.')
  image = prefix + data[0] + '.asc'
  os.system('echo scaled image = ' + image)
  os.system('./proc.scale ' + image + ' 250 1000')
  os.system("echo "+ image + " frame scaling works!")

  return

def genlat(prefix, files):

  fout = open(prefix + '/genlat.input', 'w')

  for file in files:
    os.system('cxi.image2pickle ' + prefix + '/processed/proc.' + file)
    file_prefix = file.rstrip('.img')
    scale = open(prefix + '/processed/proc.' + file_prefix + '.scale', 'r')
    lines = scale.readlines()
    scale_new = lines[0]
    scale_newest = scale_new.rstrip('\n')
    scale.close()

    fout.write(prefix + '/processed/proc.' + file_prefix + '.pickle ' + scale_newest + '\n')
  
  fout.write('.')
  fout.close()

  return

def args_generator(location_diffuse, index_1, index_2, index_3, cella, cellb, cellc, resolution, lattice_name, processors, u_p, k_setting, input):

  #Will need a way to read in the known_setting parameter

  os.system('cxi.image2pickle ' + location_diffuse + '/' + index_1)
  os.system('cxi.image2pickle ' + location_diffuse + '/' + index_2)
  os.system('cxi.image2pickle ' + location_diffuse + '/' + index_3)


  split_1 = index_1.split('.')
  location_1 = split_1[0]

  split_2 = index_2.split('.')
  location_2 = split_2[0]

  split_3 = index_3.split('.')
  location_3 = split_3[0]

  args = ''

  args += 'indexing.data=%s/%s.pickle ' %(location_diffuse, location_1)
  args += 'indexing.data=%s/%s.pickle ' %(location_diffuse, location_2)
  args += 'indexing.data=%s/%s.pickle ' %(location_diffuse, location_3)
  args += 'codecamp.maxcell=800 '
  args += 'index_only=True '
  args += 'analyze.image=45 '
  args += 'file_format=%s ' %input
  args += 'diffuse.lattice.resolution=%0.2f ' %resolution
  args += 'cell.a=%0.2f ' %float(cella)
  args += 'cell.b=%0.2f ' %float(cellb)
  args += 'cell.c=%0.2f ' %float(cellc)
  args += 'inputlist.fname=%s/genlat.input ' %location_diffuse
  args += 'diffuse.lattice.fname=%s ' %lattice_name
  args += 'target_cell=%s,%s,%s,%s,%s,%s ' %(u_p[0], u_p[1], u_p[2], u_p[3], u_p[4], u_p[5])
  args += 'known_setting=%d ' %int(k_setting)
  args += 'np=%d' %processors


  return args

def single_conversion(map_1):

  f_1 = open(map_1,'r')

  lines_1 = f_1.readlines()

  lattice = dict()

  for line in lines_1:
    signal = line.split()
    if len(signal) == 4:
      h = int(signal[0])
      k = int(signal[1])
      l = int(signal[2])
      intensity = float(signal[3])

      if h not in lattice:
         lattice[h] = dict()

      if k not in lattice[h]:
         lattice[h][k] = dict()

      if l not in lattice[h][k]:
          lattice[h][k][l] = dict()

          lattice[h][k][l]["Signal_1"] = float(intensity)
          #print intensity
  return lattice

def map_symmetry_extension(map, unitcell):
###Returns symmetrized mtz format map

  #Convert vtk map to lat and then to hkl
  os.system('vtk2lat ' + map + ' raw_output.lat')
  os.system('lat2hkl raw_output.lat raw_output.hkl')


  #Read in hkl file and populate miller array
  from cctbx import crystal
  inf = open('raw_output.hkl', "r")
  indices = flex.miller_index()
  i_obs = flex.double()
  #sig_i = flex.double()
  for line in inf.readlines():
    assert len(line.split())==4
    line = line.strip().split()
    i_obs_ = float(line[3])#/10000 #10000 is a uniform scale factor meant to re-size all diffuse intensities (normally too large for scalepack)
    #sig_i_ = math.sqrt(i_obs_) 
    #if(abs(i_obs_)>1.e-6): # perhaps you don't want zeros
    indices.append([int(line[0]),int(line[1]),int(line[2])])
    i_obs.append(i_obs_)
    #sig_i.append(sig_i_)
  inf.close()

  #Get miller array object
  ###NEED TO SET MTZ ARRAY "EXPAND TO P1"
  cs = crystal.symmetry(unit_cell=(float(args[1]), float(args[2]), float(args[3]), float(args[4]), float(args[5]), float(args[6])), space_group=args[7])
  miller_set=miller.set(cs, indices, anomalous_flag=False)
  ma = miller.array(miller_set=miller_set, data=i_obs, sigmas=None)
  ma.set_observation_type_xray_intensity()
  mtz_dataset = ma.as_mtz_dataset(column_root_label="Intensity")
  mtz_dataset.mtz_object().write('output_p1.mtz')

  return final

def friedel_hkl(map):
###Converts mtz to hkl format map and applies Friedel's law

  os.system('phenix.mtz.dump -c output_p1.mtz > output_p1.hkl')

  fin = open(map + 'output_p1.hkl','r')
  fout = open(map+ 'output_friedel.hkl', 'w')

  lines = fin.readlines()

  for line in lines[10:]:


    data = line.split()

    if len(data) == 4:

      h = int(data[0])
      h_new = -1*h
      k = int(data[1])
      k_new = -1*k
      l = int(data[2])
      l_new = -1*l
      i = float(data[3])

      #sig = data[4]

      x_1 = (str(h)+ ' ' + str(k) + ' ' + str(l) + ' ' + str(i) + '\n')
      x_2 = (str(h_new)+ ' ' + str(k_new) + ' ' + str(l_new) + ' ' + str(i) + '\n')

      fout.write(x_1)
      fout.write(x_2)

  fin.close()
  fout.close()

def hkl_to_vtk(file, cella, cellb, cellc, res, file_name):
  #Origin is defined at 0 0 0
  vtkfile = open(file_name + '.vtk', 'w')
  latxdim = (int(cella/res)+1)*2
  latydim = (int(cellb/res)+1)*2
  latzdim = (int(cellc/res)+1)*2

  a_recip = 1./cella
  b_recip = 1./cellb
  c_recip = 1./cellc

  i0=latxdim/2-1
  j0=latydim/2-1
  k0=latzdim/2-1
  latsize = latxdim*latydim*latzdim

  print >>vtkfile,"# vtk DataFile Version 2.0"
  print >>vtkfile,"Generated using labelit tools"
  print >>vtkfile,"ASCII"
  print >>vtkfile,"DATASET STRUCTURED_POINTS"
  print >>vtkfile,"DIMENSIONS %d %d %d"%(latxdim,latydim,latzdim)
  print >>vtkfile,"SPACING %f %f %f"%(a_recip,b_recip,c_recip)
  print >>vtkfile,"ORIGIN %f %f %f" %(-i0*a_recip,-j0*b_recip,-k0*c_recip)
  print >>vtkfile,"POINT_DATA %d"%(latsize)
  print >>vtkfile,"SCALARS volume_scalars float 1"
  print >>vtkfile,"LOOKUP_TABLE default\n"

    #print data[1][1][1]["Signal_1"]

  for i in range(-num_a,num_a+1):
    for j in range(-num_b,num_b+1):
      for k in range(-num_c,num_c+1):
        if i in data and j in data[i] and k in data[i][j]:      
          print >>vtkfile,data[i][j][k]["Signal_1"], 
          #print 'hello!!!'
        else:
          print >>vtkfile,"0.0", 


      print >>vtkfile,""


def aniso_convert(file_name, cell_a, cell_b, cell_c, resolution, file_format):

  os.system('hkl2vtk output_p1.hkl output_p1.vtk template.vtk')
  os.system('vtk2lat output_p1.vtk output.lat')
  os.system('avgrlt output.lat output.rf')
  os.system('subrflt output.rf output.lat anisotropic.lat')
  os.system('lat2hkl anisotropic.lat anisotropic.hkl')

  #Read in hkl file and populate miller array
  from cctbx import crystal
  inf = open('anisotropic.hkl', "r")
  indices = flex.miller_index()
  i_obs = flex.double()
  #sig_i = flex.double()
  for line in inf.readlines():
    assert len(line.split())==4
    line = line.strip().split()
    i_obs_ = float(line[3])#/10000 #10000 is a uniform scale factor meant to re-size all diffuse intensities (normally too large for scalepack)
    #sig_i_ = math.sqrt(i_obs_) 
    #if(abs(i_obs_)>1.e-6): # perhaps you don't want zeros
    indices.append([int(line[0]),int(line[1]),int(line[2])])
    i_obs.append(i_obs_)
    #sig_i.append(sig_i_)
  inf.close()

  #Get miller array object
  cs = crystal.symmetry(unit_cell=(float(args[1]), float(args[2]), float(args[3]), float(args[4]), float(args[5]), float(args[6])), space_group=args[7])
  miller_set=miller.set(cs, indices, anomalous_flag=False)
  ma = miller.array(miller_set=miller_set, data=i_obs, sigmas=None)
  ma.set_observation_type_xray_intensity()
  mtz_dataset = ma.as_mtz_dataset(column_root_label="Intensity")
  mtz_dataset.mtz_object().write('anisotropic.mtz')


def hkl2vtk(data, cell_a, cell_b, cell_c, resolution, name):
  #Origin is defined at 0 0 0
  vtkfile = open(name, 'w')

  num_a = int(cell_a/resolution)
  num_b = int(cell_b/resolution)
  num_c = int(cell_c/resolution)

  latsize = (2*num_a+1)*(2*num_b+1)*(2*num_c+1)

  print >>vtkfile,"# vtk DataFile Version 2.0"
  print >>vtkfile,"Generated using labelit tools"
  print >>vtkfile,"ASCII"
  print >>vtkfile,"DATASET STRUCTURED_POINTS"
  print >>vtkfile,"DIMENSIONS %d %d %d"%(num_a*2+1,num_b*2+1,num_c*2+1)
  print >>vtkfile,"SPACING %f %f %f"%(1/cell_a,1/cell_b,1/cell_c)
  print >>vtkfile,"ORIGIN %.8f %.8f %.8f" %(-num_a/cell_a,-num_b/cell_b,-num_c/cell_c)
  print >>vtkfile,"POINT_DATA %d"%(latsize)
  print >>vtkfile,"SCALARS volume_scalars float 1"
  print >>vtkfile,"LOOKUP_TABLE default\n"

    #print data[1][1][1]["Signal_1"]

  for i in range(-num_a,num_a+1):
    for j in range(-num_b,num_b+1):
      for k in range(-num_c,num_c+1):
        if i in data and j in data[i] and k in data[i][j]:      
          print >>vtkfile,data[i][j][k]["Signal_1"], 
          #print 'hello!!!'
        else:
          print >>vtkfile,"0.0", 


      print >>vtkfile,""



def vtk2hkl(file, cell_a, cell_b, cell_c, resolution):

  fin = open(file,'r')
  fout = open('aniso.hkl','w')

  lines = fin.readlines()

  num_a = int(cell_a/resolution)
  num_b = int(cell_b/resolution)
  num_c = int(cell_c/resolution)

  h = -num_a
  k = -num_b
  print h
  print k

  count = 1

  for line in lines[10:]:
    l = -num_c
    #print l
    data = line.split()
    print data
    dic = dict()

    count_item=0


    for item in data:
      if h not in dic:
        dic[h] = dict()
      if k not in dic[h]:
        dic[h][k] = dict()
      if l not in dic[h][k]:
        dic[h][k][l] = dict()
      
      dic[h][k][l]["Signal"] = float(item)
      #print h,k,l, float(item)
      l+=1
      
    k += 1

    if k == num_b+1:
      k = -num_b
      h += 1


  for key_h in dic:
    print key_h
    for key_k in dic[key_h]:
      for key_l in dic[key_h][key_k]:
        fout.write(str(key_h) + ' ' + str(key_k) + ' ' + str(key_l) + ' ' + str(dic[key_h][key_k][key_l]["Signal"]) + '\n')

if __name__=="__main__":

  import sys
  import os
  import math
  
  run(sys.argv[1])


