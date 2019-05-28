## RSQlite Practices

mydb <- dbConnect(RSQLite::SQLite(), "my-db.sqlite")
dbDisconnect(mydb)
unlink("my-db.sqlite")


