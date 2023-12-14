import pandas as pd

df = pd.read_csv('foldseek_result_cluster.tsv', sep = '\t')

df.columns = ['Representative', 'Clustered']

value_counts = df['Representative'].value_counts().to_dict()

count_keys_with_value_one = sum(1 for key, value in value_counts.items() if value == 1)

print(count_keys_with_value_one)

## 68 singletons in the cluster, not sure whether worth it to drop or not, may cause unmapped values when analyzing ##
