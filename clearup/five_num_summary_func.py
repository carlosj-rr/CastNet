def five_num_summary(dataset):
     mini=np.min(dataset)
     maxi=np.max(dataset)
     dset_len=len(dataset)
     sdataset=sorted(dataset)
     if dset_len%2 != 0:
             median_idx=int((dset_len+1)/2) -1
             data_for_q1 = sdataset[:median_idx]
             data_for_q3 = sdataset[median_idx:]
     else:   
             data_for_q1 = sdataset[:int((dset_len/2))]
             data_for_q3 = sdataset[int(dset_len/2):]
     q1 = np.median(data_for_q1)
     q3 = np.median(data_for_q3)
     median=np.median(sdataset)
     iqr=q3-q1
     return {'min': mini, 'Q1': q1, 'Q2(median)': median, 'Q3': q3, 'Q4(max)': maxi, 'IQR': iqr}
