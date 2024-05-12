###Filtering the human RNA-binding site dataset provided by POSTAR3 to exclusively contain DDX6 binding sites###

import sqlite3
import pandas as pd
import csv

#To start with! The file 'human_RBP_binding_sites.txt' is ~4.5gb- fairly big.
#There are no headers, but there should be 10 columns.
#Values are separated by '\t'
#if using pandas, may need to separate into chunks, as the below process was killed when trying to process entire dataset. First: try using sqlite

#df = pd.read_csv('human_RBP_binding_sites.txt', header = None, sep = '\t')
#print(df.head())

#for sqlite: initially create a database and a table. We'll need to transfer the data from the txt file to the table.

conn = sqlite3.connect('RBP_binding_sites.db')

#create a cursor, which is what tells the database what you want to do
c = conn.cursor()

#create a table- DONE
#c.execute("""CREATE TABLE RBPs (
 #       chromosome TEXT,
 #       start_coord INTEGER,
 #       end_coord INTEGER,
 #       peak_id TEXT,
 #       strand TEXT,
 #       RBP_name TEXT,
 #       experiment_method TEXT,
 #       sample_or_tissue_used TEXT,
 #       accession_of_raw_data TEXT,
 #       confidence_score REAL
 #   )""")

#table created; now need to transfer data into table- DONE

#with open('human_RBP_binding_sites.txt', 'r') as fin:
#    reader = csv.reader(fin, delimiter='\t')  # Specify space as the delimiter
#    for row in reader:
 #       # Assuming the order of values matches the order of placeholders in your SQL statement
#        if len(row) == 10:  # Ensure that the row has exactly 10 values
#            chromosome, start_coord, end_coord, peak_id, strand, RBP_name, experiment_method, sample_or_tissue_used, accession_of_raw_data, confidence_score = row
#            conn.execute("INSERT INTO RBPs (chromosome, start_coord, end_coord, peak_id, strand, RBP_name, experiment_method, sample_or_tissue_used, accession_of_raw_data, confidence_score) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);",
#                         (chromosome, start_coord, end_coord, peak_id, strand, RBP_name, experiment_method, sample_or_tissue_used, accession_of_raw_data, confidence_score))
#        else:
#            print(f"Skipping row with incorrect number of values: {row}")

#THE ABOVE WORKED! SQLITE TABLE 'RBPs' SHOULD NOW CONTAIN THE TEXT FILE VALUES

#Query the database
c.execute("SELECT * FROM RBPs")
print(c.fetchone())
#IT WORKS

#now want to specifically query for all rows containing 'DDX6' under 'RBP_name'
RBP = 'DDX6'
c.execute("""
    SELECT t1.*
    FROM RBPs t1
    INNER JOIN (SELECT '%' || ? || '%' AS RBP) t2
    ON t1.RBP_name LIKE t2.RBP
""", (RBP,))

# Fetch the filtered rows
filtered_rows = c.fetchall()

# Process the filtered data as needed
#for row in filtered_rows:
#    print(row)

#The above confirms that the query worked. Now need to export the file as a .csv.
# Specify the output CSV file path
csv_file_path = 'DDX6_binding_coords.csv'

# Write the filtered data to the CSV file- DONE
#with open(csv_file_path, 'w', newline='') as csvfile:
#    csv_writer = csv.writer(csvfile)
    # Write header (column names)
#    csv_writer.writerow([col[0] for col in c.description])
    # Write data rows
#    csv_writer.writerows(filtered_rows)

#print(f"Filtered data exported to {csv_file_path}")


#commit our command
conn.commit()

#close our connection
conn.close()
