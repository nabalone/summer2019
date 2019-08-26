import pandas as pd
db = pd.read_table('alertstable_v3',sep=None,index_col = False,
			   engine='python')
db2 = pd.read_table('alertstable_v3.lasthalf',sep=None,index_col = False,
				engine='python')
db = db.append(db2,ignore_index=True)

idNum = 60270
print(db.where(db['eventID'] == idNum).dropna())
	
	