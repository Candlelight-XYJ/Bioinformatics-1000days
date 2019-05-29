## RSQlite Practices
library(RSQLite)

mydb <- dbConnect(RSQLite::SQLite(), "my-db.sqlite")
dbDisconnect(mydb)
unlink("my-db.sqlite")

## create temporary database using ""
mydb <- dbConnect(RSQLite::SQLite(), "")
dbDisconnect(mydb)

## Loading data
mydb <- dbConnect(RSQLite::SQLite(), "")
dbWriteTable(mydb, "mtcars", mtcars)
dbWriteTable(mydb, "iris", iris)
dbListTables(mydb)
#> [1] "iris"   "mtcars"

## Queries
dbGetQuery(mydb, 'SELECT * FROM mtcars LIMIT 5')
dbGetQuery(mydb, 'SELECT * FROM iris WHERE "Sepal.Length" < 4.6')

dbGetQuery(mydb, 'SELECT * FROM iris WHERE "Sepal.Length" < :x', 
           params = list(x = 4.6))

## Batched queries
rs <- dbSendQuery(mydb, 'SELECT * FROM mtcars')
while (!dbHasCompleted(rs)) {
  df <- dbFetch(rs, n = 10)
  print(nrow(df))
}
#> [1] 10
#> [1] 10
#> [1] 10
#> [1] 2
dbClearResult(rs)


## Multiple parameterised queries

# Call dbBind() to set the parameters:
rs <- dbSendQuery(mydb, 'SELECT * FROM iris WHERE "Sepal.Length" < :x')
dbBind(rs, param = list(x = 4.5))
nrow(dbFetch(rs))
#> [1] 4
dbBind(rs, param = list(x = 4))
nrow(dbFetch(rs))
#> [1] 0
dbClearResult(rs)


# pass multiple parameters in one call to dbBind():
rs <- dbSendQuery(mydb, 'SELECT * FROM iris WHERE "Sepal.Length" = :x')
dbBind(rs, param = list(x = seq(4, 4.4, by = 0.1)))
nrow(dbFetch(rs))
#> [1] 4
dbClearResult(rs)


## statement

dbExecute(mydb, 'DELETE FROM iris WHERE "Sepal.Length" < 4')
#> [1] 0
rs <- dbSendStatement(mydb, 'DELETE FROM iris WHERE "Sepal.Length" < :x')
dbBind(rs, param = list(x = 4.5))
dbGetRowsAffected(rs)
#> [1] 4
dbClearResult(rs)


## statement continue
con <- dbConnect(RSQLite::SQLite(), ":memory:")
dbWriteTable(con, "cars", head(cars, 3))
rs <- dbSendStatement(
  con,
  "INSERT INTO cars (speed, dist) VALUES (1, 1), (2, 2), (3, 3)"
)
dbHasCompleted(rs)
dbGetRowsAffected(rs)
dbClearResult(rs)
dbReadTable(con, "cars")   # there are now 6 rows

# Pass one set of values directly using the param argument:
rs <- dbSendStatement(
  con,
  "INSERT INTO cars (speed, dist) VALUES (?, ?)",
  param = list(4L, 5L)
)
dbClearResult(rs)

# Pass multiple sets of values using dbBind():
rs <- dbSendStatement(
  con,
  "INSERT INTO cars (speed, dist) VALUES (?, ?)"
)
dbBind(rs, list(5:6, 6:7))
dbBind(rs, list(7L, 8L))
dbClearResult(rs)
dbReadTable(con, "cars")   # there are now 10 rows

dbDisconnect(con)

