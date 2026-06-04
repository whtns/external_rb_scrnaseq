import numpy as np
import pandas as pd
from pandas import DataFrame, Series
import sqlite3 as db

type_diffex = pd.read_table("/dataVolume/storage/scEiad/human_pseudobulk/diff_resultsCellType.tsv.gz", engine = "pyarrow")

database = "/dataVolume/storage/scEiad/human_pseudobulk/diff_resultsCellType.sqlite"
conn = db.connect(database)

type_diffex.to_sql(name='diffex', con=conn)
conn.close()

