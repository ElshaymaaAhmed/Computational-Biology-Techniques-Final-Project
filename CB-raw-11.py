#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pyopenms
from pyopenms import *


# In[2]:


exp = MSExperiment()
MzMLFile().load("D:/projects/Computional Biology/Fusion_180220_11.mzML",exp)
spectra = exp.getSpectra()


# In[12]:


observed_spectrum = spectra[6]
print(observed_spectrum)


# ## Proteolytic Digestion with Trypsin

# In[13]:


dig = ProteaseDigestion()
dig.getEnzymeName() 
bsa = "".join([l.strip() for l in open("Scerevisiae_UPS2_1802.fasta").readlines() if l.startswith('>') == False])
bsa = AASequence.fromString(bsa)
result = []
dig.digest(bsa, result)
len(result) 


# In[14]:


print(result[6].toString())


# ## Now we generate the theoretical spectrum of that peptide: 

# In[19]:


tsg = TheoreticalSpectrumGenerator()
theo_spectrum = MSSpectrum()
p = tsg.getParameters()
p.setValue("add_y_ions", "true")
p.setValue("add_b_ions", "true")              
p.setValue("add_metainfo", "true")                          
tsg.setParameters(p)   

peptide = result[6]
tsg.getSpectrum(theo_spectrum, peptide, 1, 2)


# ## Now we can plot the observed and theoretical spectrum as a mirror plot:
# 

# In[20]:


import numpy as np
from matplotlib import pyplot as plt

def mirror_plot(obs_mz, obs_int, theo_mz, theo_int, title):
  
    obs_int = [element / max(obs_int) for element in obs_int]
    theo_int = [element * -1 for element in theo_int] 
    plt.figure(figsize=(12,8))
    plt.bar(obs_mz, obs_int, width = 3.0)
    plt.bar(theo_mz, theo_int, width = 3.0)
    plt.title(title)
    plt.ylabel('intensity')
    plt.xlabel('m/z')

obs_mz, obs_int = observed_spectrum.get_peaks()

print(min(obs_mz))
print(max(obs_mz))

theo_mz, theo_int = [], []
for mz, intensity in zip(*theo_spectrum.get_peaks()):
    if mz >= 200.0 and mz <= 800.0:
        theo_mz.append(mz)
        theo_int.append(intensity)

title = 'Observed vs theoretical spectrum'
mirror_plot(obs_mz, obs_int, theo_mz, theo_int, title)


# ## Now we want to find matching peaks between observed and theoretical spectrum. 

# In[21]:


alignment = []
spa = SpectrumAlignment()
p = spa.getParameters()
p.setValue("tolerance", 0.5)
p.setValue("is_relative_tolerance", "false")
spa.setParameters(p)
spa.getSpectrumAlignment(alignment, theo_spectrum, observed_spectrum)


# ## The alignment contains a list of matched peak indices. We can simply inspect matching peaks with: 

# In[22]:


print("Number of matched peaks: " + str(len(alignment)))
print("ion\ttheo. m/z\tobserved m/z")

for theo_idx, obs_idx in alignment:
    ion_name = theo_spectrum.getStringDataArrays()[0][theo_idx].decode()
    ion_charge = theo_spectrum.getIntegerDataArrays()[0][theo_idx]
    print(ion_name + "\t" + str(ion_charge) + "\t"
      + str(theo_spectrum[theo_idx].getMZ())
      + "\t" + str(observed_spectrum[obs_idx].getMZ()))


# ## The mirror plot can also be used to visualize the aligned spectrum: 

# In[23]:


theo_mz, theo_int, obs_mz, obs_int = [], [], [], []
for theo_idx, obs_idx in alignment:
    theo_mz.append(theo_spectrum[theo_idx].getMZ())
    theo_int.append(theo_spectrum[theo_idx].getIntensity())
    obs_mz.append(observed_spectrum[obs_idx].getMZ())
    obs_int.append(observed_spectrum[obs_idx].getIntensity())

title = 'Observed vs theoretical spectrum (aligned)'
mirror_plot(obs_mz, obs_int, theo_mz, theo_int, title)


# In[ ]:




