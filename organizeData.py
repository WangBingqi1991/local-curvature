import scipy.io
import numpy as np
import pandas as pd


mat_data = scipy.io.loadmat('C:/Users/DAG/Desktop/bro/1025data/finaml1025image_all_stats.mat')
allStats = mat_data['allStats']

data_rows = []

for j in range(allStats.shape[0]):
    for i in range(allStats.shape[1]):
        cell_data = allStats[j, i]
        d9595 = cell_data['d9595'][0][0][0][0]
        print(d9595)
        numValidVesicles = cell_data['numValidVesicles'][0][0][0][0]

        if 'vesicleDetails' in cell_data.dtype.names:
            vesicle_details = cell_data['vesicleDetails'][0, 0]

            print(vesicle_details.size)
            for k in range(vesicle_details.size):
                vesicle = vesicle_details[k][0]
                isValid = vesicle['isValid'][0][0]

                areaNm2 = vesicle['areaNm2'][0][0]

                perimeters = vesicle['perimeter'][0][0]

                hotspot_counts = vesicle['hotspotCounts'].flatten() if vesicle_details['hotspotCounts'].size > 0 else np.full(                    29, 0)
                hotspot_counts_padded = np.pad(hotspot_counts, (0, max(0, 29 - len(hotspot_counts))),
                                               constant_values=0)[:29]

                row = [i,k,d9595,numValidVesicles,isValid, areaNm2, perimeters] + hotspot_counts_padded.tolist()
                data_rows.append(row)

columns = ['k', 'i','d9595', 'numValidVesicles', 'isValid' , 'areaNm2', 'perimeters'] + [f'hotspotCount_{k + 1}' for k in range(29)]

df = pd.DataFrame(data_rows, columns=columns)

df = df.fillna(0)

df.to_excel('C:/Users/DAG/Desktop/bro/1025data/finaml1025image_all_stats.xlsx', index=False)