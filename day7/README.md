学习来源：
https://cran.r-project.org/web/packages/RSQLite/vignettes/RSQLite.html

[toc]

## 1. 创建一个新的SQLite数据库
我们可以通过函数 **`dbConnect()`** 来连接和创建一个新的数据库
```r
mydb <- dbConnect(RSQLite::SQLite(), "my-db.sqlite")
dbDisconnect(mydb)
unlink("my-db.sqlite")
```
可以使用双引号**`""`**来构建一个桌面的临时数据库，或者使用 **`:memory:`**和 **`file::memory`**
构建一个内存里的数据库。当我们断开连接，临时数据库就会直接自动删除。

```r
mydb <- dbConnect(RSQLite::SQLite(), "")
dbDisconnect(mydb)
```

## 2. Loading data (加载数据)
我们可以使用函数 **`dbWriteTable():`** 来将R数据框格式的数据，导入数据框
```r
mydb <- dbConnect(RSQLite::SQLite(), "")
dbWriteTable(mydb, "mtcars", mtcars)
dbWriteTable(mydb, "iris", iris)
dbListTables(mydb)
#> [1] "iris"   "mtcars"
```
## 3. Queries (查询) 
使用 **`dbGetQuery()`** 函数执行查询
```r
dbGetQuery(mydb, 'SELECT * FROM mtcars LIMIT 5')
#>    mpg cyl disp  hp drat    wt  qsec vs am gear carb
#> 1 21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
#> 2 21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
#> 3 22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
#> 4 21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
#> 5 18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
```
**`注意：`**  不要使用`paste()`函数来构建用户查询，这样很可能会引起SQL注入攻击
要使用params参数来查询
```r
dbGetQuery(mydb, 'SELECT * FROM iris WHERE "Sepal.Length" < :x', 
  params = list(x = 4.6))
#>   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> 1          4.4         2.9          1.4         0.2  setosa
#> 2          4.3         3.0          1.1         0.1  setosa
#> 3          4.4         3.0          1.3         0.2  setosa
#> 4          4.5         2.3          1.3         0.3  setosa
#> 5          4.4         3.2          1.3         0.2  setosa
```
这会显著提高数据库的安全性

## 4. Batched queries
如果我们执行了一个查询，但是查询结果很大，超出了内存范围，不能直接储存于内存中，
此时我们就需要使用 **`dbSendQuery()`** ，**`dbFetch()`** ， **`dbClearResults()`** 函数来批量处理检索结果，从超大的查询结果中选取自己想要的信息后，再确定存储下来。 
+ **`注：`** 这属于检索大文件的操作

**`dbFetch()`** 默认检索所有可能的行，使用`n`参数，可以设置返回的最大行数
```r
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
```
+ **`dbSendQuery() `** 可以执行一条查询，但是它不会直接抓取数据库中的数据记录
+ 使用 **`dbFetch() `** 对 **`dbSendQuery() `** 查询得到的object，进行数据的抓取
+ 完成 **`dbSendQuery()`**  和 **`dbFetch()`**操作后，使用 **`dbClearResult()`** 来停止抓取数据的操作。

##  5. Multiple parameterised queries

**`设置一个参数，进行参数化查询`**
我们也可以使用 **`dbBind()`** 函数，结合 **`dbSendQuery()`** 函数以及**`dbFetch()`** 函数来完成参数化查询
```r
rs <- dbSendQuery(mydb, 'SELECT * FROM iris WHERE "Sepal.Length" < :x')
dbBind(rs, param = list(x = 4.5))
nrow(dbFetch(rs))
#> [1] 4
dbBind(rs, param = list(x = 4))
nrow(dbFetch(rs))
#> [1] 0
dbClearResult(rs)
```

**`设置多个参数，进行参数化查询`**
```r
rs <- dbSendQuery(mydb, 'SELECT * FROM iris WHERE "Sepal.Length" = :x')
dbBind(rs, param = list(x = seq(4, 4.4, by = 0.1)))
nrow(dbFetch(rs))
#> [1] 4
dbClearResult(rs)
```

## 6. Satements
DBI有两个新的函数dbSendStatement()和dbExecute()，它们是dbSendQuery()和dbGetQuery()的对应函数，应用于不返回表结果的SQL语句查询（例如将记录插入表、更新表或设置引擎参数）。如果不是很需要得到返回的查询结果，使用新函数是一种很好的实践。
```r
dbExecute(mydb, 'DELETE FROM iris WHERE "Sepal.Length" < 4')
#> [1] 0
rs <- dbSendStatement(mydb, 'DELETE FROM iris WHERE "Sepal.Length" < :x')
dbBind(rs, param = list(x = 4.5))
dbGetRowsAffected(rs)
#> [1] 4
dbClearResult(rs)
```

```r
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
```

## 参考连接
+ https://www.runoob.com/sqlite/sqlite-drop-table.html
+ https://www.runoob.com/sqlite/sqlite-drop-table.html
+ https://cran.r-project.org/web/packages/RSQLite/vignettes/RSQLite.html