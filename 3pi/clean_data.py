import pandas as pd

filenames_to_clean = ['sn1a.csv', 'sn1b.csv', 'sn1c.csv', 'snii.csv', 'sniin.csv']
types = [('Ia',), ('Ib', 'Ib/c'), ('Ic', 'Ib/c'), ('II',), ('IIn',)]
def clean(filename, type):
        data = pd.read_csv(filename)
        data = data.loc[:, ['Name', 'Type', 'R.A.', 'Dec.', 'Phot.']]
        data = data.where(data['Type'] == type).dropna()
        data = data.where(data['Phot.'] > 30).dropna()
        data.to_csv(filename[:-4] + '_clean.csv')