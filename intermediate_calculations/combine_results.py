#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import glob


# In[2]:


input_list = glob.glob()
cc_thresh = 
output_file = 


# In[3]:


frame = pd.read_csv(input_list[0])

file_list = [None]*len(input_list)
count = 1

for i in range(len(input_list)):
    
    frame = pd.read_csv(input_list[i])
    
    frame = frame[frame['cc'] > cc_thresh]
    
    file_list[i] = frame
    #print(frame.shape[0], str(count)+' of '+ str(len(input_list)))
    
    count += 1


# In[4]:


df = pd.concat(file_list, axis = 0, ignore_index = True)

aggregation_functions = {'cc':'mean', 'AR':'mean', 'cha':list}
df_new = df.groupby(['dtm', 'dtn'], as_index=False).aggregate(aggregation_functions).reindex(columns=df.columns)


# In[5]:


df_new.to_csv(output_file)


# In[ ]:




