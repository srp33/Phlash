import sqlite3

def merge0_0(DATABASE):
    # Create a connection object
    connection  = sqlite3.connect(DATABASE)

    # Get a cursor
    cursor      = connection.cursor()

    # Make necessary changes to tables
    dropTable = 'DROP TABLE metadata'
    createTable = 'CREATE TABLE metadata (version INTEGER, PRIMARY KEY(version));'
    addColumn = 'ALTER TABLE settings ADD COLUMN prodigal BOOLEAN DEFAULT "True" NOT NULL'

    # # Execute changes
    try:
        cursor.execute(dropTable)
    finally:
        cursor.execute(createTable)
    cursor.execute(addColumn)

    # Query the SQLite master table
    tableQuery = "select * from sqlite_master"
    cursor.execute(tableQuery)
    tableList = cursor.fetchall()

    # Print the updated listed of tables after renaming the stud table
    for table in tableList:
        print("Database Object Type: %s"%(table[0]))
        print("Name of the database object: %s"%(table[1]))
        print("Name of the table: %s"%(table[2]))
        print("Root page: %s"%(table[3]))
        print("SQL Statement: %s"%(table[4]))

    # close the database connection
    connection.close()

    return "Merged Successfully"