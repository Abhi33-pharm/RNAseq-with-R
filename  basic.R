#String or character
#Find the Length of a character

strg<- "TGCA"
nchar(strg)

# Join two strings 

s1<- "Naem"
s2<-"Islam"

str<- paste(s1,s2)
str
nchar(str)

#Compare the string 

s1<- "TCGA"
s2<- "TCGA"

# Changing string case 

tolower("TCGA")

nchar(s1)
nchar(s2)

s1==s2



# Creating Vector Using C() function 

num<- c(1,3,4)
str<- c("black","red", "blue")
is_holiday<- c(T, F, T,F)
mixed<- c(2,"black", T)


class(num)
class(str)
class(is_holiday)
class(mixed)


# Items or Elements
# Position 

str<- c("black","red", "blue")
str[2]

# Create Vector Using ':' operator 
seq_vector<- 1:5
seq_vector

# Create Vector Using seq() function
# seq(start=, end=, step=)

seq(1,15)
seq(1,15,3)
ages<- c(12,35,45,78,34,67)

#Sub-setting elements 

ages[3]
ages[2:6]
ages[c(2,5,6)]

# Factors (levels= number of groups, obs)

is_smoking<- factor(c("yes","no","no","no","yes"))
is_smoking

# Matrix
matrix(1:9)
matrix(1:9,ncol = 3)
matrix(1:9,ncol = 3,nrow = 3)
mat<-matrix(1:9,ncol = 3,nrow = 3,byrow = T)

#matrix properties
dim(mat)
ncol(mat)
nrow(mat)

# sub-setting

mat[1]
mat[2,3]

#Accessing entire row and column 

mat[2,] #row
mat[,2] #column 


#List 

my_list<- list(
  name= "abhi",
  vec= c(2,4,5),
  mat=matrix(1:9,ncol = 3,nrow = 3,byrow = T),
  is_smoking=factor(c("yes","no","no","no","yes"))
)
my_list


#Accessing the list 
my_list$is_smoking


# Exploring Build-in data
data()
Titanic

#Data frame

df<- data.frame(
  age= c(22,34,78,30),
  sex= c("M","F","M","F"),
  smoking_status= c("yes","no","no","no")
)
df

#Exploring data of the data frame

df$age
df$sex
df$smoking_status

#Checking data type
str(df)

#Conversion of the data type (as.family function)
df$sex<-as.factor(df$sex)
df$smoking_status<-as.factor(df$smoking_status)

# If-else statement: If (condition) {make statement} : for 2 conditions

num<- -10

if (num>0){
  print("positive")
}else {
  print("negative")
}

#if...else if...else...if...else (more than 2 conditions)
if (condition){
  #make a decision
} else if (condition){
  #make a decision
} else if (condition){
  #make a decision
} else {
  
}

bmi<- 30

if (bmi<18.5){
  print("underweight")
} else if (bmi>=18.5 & bmi<25){
  print("normal")
} else if (bmi>= 25 & bmi<30){
  print("overweight")
}else {
  print("obese")
}

# By using if else function

num<- 20
ifelse(num>0, "positive", "negative")

ages<- c(12,35,45,78,34,67)

ifelse(ages<30, "young","older")
