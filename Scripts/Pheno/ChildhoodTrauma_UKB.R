##### Childhood trauma definition in line with GLAD #####
##### Rujia Wang 20241104 ####
library(dplyr)
library(ukbkings)
library(dplyr)
library(stringr)

project_dir <- "/datasets/ukbiobank/ukb82087"

f <- bio_field(project_dir)

head(f)
glimpse(f)

f %>%
  distinct(basket)

f %>%
  select(field, name) %>%
  filter(str_detect(field, "20487|20488|20489|20490|20491|2052|20530|20531|29076|29077|29078|29079|2908|29090")) %>%
  bio_field_add("CT_MHQ.txt")

bio_phen(
  project_dir,
  field = "CT_MHQ.txt",
  out = "CT_MHQ"
)

trauma_data<-readRDS("CT_MHQ.rds")

#### MHQ1 ####

#### Emotional abuse: 20487 - when I was growing up, I felt that someone in my family hated me ####
#### 0= Never true, 1= Rarely true, 2= Sometimes true, 3=Often, 4=Very often true ####

trauma_data <- trauma_data %>%
  mutate(
    EA_MHQ1 = case_when(
      `20487-0.0` > 1 ~ 1,             # Having emotional abuse if the answer is sometimes true, often and very often true
      `20487-0.0` == 0 | `20487-0.0` == 1 ~ 0, # No emotional abuse if the value is 0 or 1
      `20487-0.0` == -818 ~ NA_real_   # Assign NA for the value -818
    )
  )

#### Physical abuse: 20488 - when I was growing up, people in my family hit me so hard that it left me with bruises or marks ####
#### 0= Never true, 1= Rarely true, 2= Sometimes true, 3=Often, 4=Very often true ####

trauma_data <- trauma_data %>%
  mutate(
    PA_MHQ1 = case_when(
      `20488-0.0` > 1 ~ 1,             # Having physical abuse if the answer is sometimes true, often and very often true
      `20488-0.0` == 0 | `20488-0.0` == 1 ~ 0, # No physical abuse if the value is 0 or 1
      `20488-0.0` == -818 ~ NA_real_   # Assign NA for the value -818
    )
  )

#### Sexual abuse: 20490 - when I was growing up, someone molested me (sexually) ####
#### 0= Never true, 1= Rarely true, 2= Sometimes true, 3=Often, 4=Very often true ####

trauma_data <- trauma_data %>%
  mutate(
    SA_MHQ1 = case_when(
      `20490-0.0` > 0 ~ 1,             # Having sexual abuse if the answer is rarely true, sometimes true, often and very often true
      `20490-0.0` == 0 ~ 0, # No sexual abuse if the value is never true
      `20490-0.0` == -818 ~ NA_real_   # Assign NA for the value -818
    )
  )

#### Emotional neglect: 20489 - when I was growing up, I felt loved ####
#### 0= Never true, 1= Rarely true, 2= Sometimes true, 3=Often, 4=Very often true ####

trauma_data <- trauma_data %>%
  mutate(
    EN_MHQ1 = case_when(
      `20489-0.0` == 0 | `20489-0.0` == 1 ~ 1, # Having emotional neglect if the answer is never true or rarely true
      `20489-0.0` > 1 ~ 0,               # No emotional neglect if the value is sometimes true, often and very often true
      `20489-0.0` == -818 ~ NA_real_   # Assign NA for the value -818
    )
  )

#### Physical neglect: 20491 - when I was growing up, there was someone to take me to the doctor if I needed it ####
#### 0= Never true, 1= Rarely true, 2= Sometimes true, 3=Often, 4=Very often true ####

trauma_data <- trauma_data %>%
  mutate(
    PN_MHQ1 = case_when(
      `20491-0.0` == 0 | `20491-0.0` == 1 ~ 1, # Having physical neglect if the answer is never true or rarely true
      `20491-0.0` > 1 ~ 0,               # No physical neglect if the value is sometimes true, often and very often true
      `20491-0.0` == -818 ~ NA_real_   # Assign NA for the value -818
    )
  )
  
#### Childhood trauma: having at least one subtype of childhood trauma defining having childhood trauma ####

trauma_data <- trauma_data %>%
  mutate(
    ChT_MHQ1 = case_when(
      PN_MHQ1 == 1 | EN_MHQ1 == 1 | EA_MHQ1 == 1 | PA_MHQ1 == 1 | SA_MHQ1 == 1 ~ 1, # Assign 1 if any of the specified columns are 1
      is.na(ChT_MHQ1) & 
      (PN_MHQ1 == 0 | EN_MHQ1 == 0 | EA_MHQ1 == 0 | PA_MHQ1 == 0 | SA_MHQ1 == 0) ~ 0 # Assign 0 if ChT_MHQ1_test is NA and any specified columns are 0
    )
  )

#### MHQ2 ####

#### Emotional abuse: 29078 - when I was growing up, I felt that someone in my family hated me ####
#### 0= Never true, 1= Rarely true, 2= Sometimes true, 3=Often true, 4=Very often true ####

trauma_data <- trauma_data %>%
  mutate(
    EA_MHQ2 = case_when(
      `29078-0.0` > 1 ~ 1,             # Having emotional abuse if the answer is sometimes true, often and very often true
      `29078-0.0` == 0 | `29078-0.0` == 1 ~ 0, # No emotional abuse if the value is 0 or 1
      `29078-0.0` == -3 ~ NA_real_   # Assign NA for the value -3
    )
  )

#### Physical abuse: 29077 - when I was growing up, people in my family hit me so hard that it left me with bruises or marks ####
#### 0= Never true, 1= Rarely true, 2= Sometimes true, 3=Often true, 4=Very often true ####

trauma_data <- trauma_data %>%
  mutate(
    PA_MHQ2 = case_when(
      `29077-0.0` > 1 ~ 1,             # Having physical abuse if the answer is sometimes true, often and very often true
      `29077-0.0` == 0 | `29077-0.0` == 1 ~ 0, # No physical abuse if the value is 0 or 1
      `29077-0.0` == -3 ~ NA_real_   # Assign NA for the value -3
    )
  )

#### Sexual abuse: 29079 - when I was growing up, someone molested me (sexually) ####
#### 0= Never true, 1= Rarely true, 2= Sometimes true, 3=Often true, 4=Very often true ####

trauma_data <- trauma_data %>%
  mutate(
    SA_MHQ2 = case_when(
      `29079-0.0` > 0 ~ 1,             # Having sexual abuse if the answer is rarely true, sometimes true, often true and very often true
      `29079-0.0` == 0 ~ 0, # No sexual abuse if the value is never true
      `29079-0.0` == -3 ~ NA_real_   # Assign NA for the value -3
    )
  )

#### Emotional neglect: 29076 - when I was growing up, I felt loved ####
#### 0= Never true, 1= Rarely true, 2= Sometimes true, 3=Often true, 4=Very often true ####

trauma_data <- trauma_data %>%
  mutate(
    EN_MHQ2 = case_when(
      `29076-0.0` == 0 | `29076-0.0` == 1 ~ 1, # Having emotional neglect if the answer is never true or rarely true
      `29076-0.0` > 1 ~ 0,               # No emotional neglect if the value is sometimes true, often and very often true
      `29076-0.0` == -3 ~ NA_real_   # Assign NA for the value -3
    )
  )

#### Physical neglect: 29080 - when I was growing up, there was someone to take me to the doctor if I needed it ####
#### 0= Never true, 1= Rarely true, 2= Sometimes true, 3=Often true, 4=Very often true ####

trauma_data <- trauma_data %>%
  mutate(
    PN_MHQ2 = case_when(
      `29080-0.0` == 0 | `29080-0.0` == 1 ~ 1, # Having physical neglect if the answer is never true or rarely true
      `29080-0.0` > 1 ~ 0,               # No physical neglect if the value is sometimes true, often and very often true
      `29080-0.0` == -3 ~ NA_real_   # Assign NA for the value -3
    )
  )
  
#### Childhood trauma: having at least one subtype of childhood trauma defining having childhood trauma ####

trauma_data$ChT_MHQ2<-NA

trauma_data <- trauma_data %>%
  mutate(
    ChT_MHQ2 = case_when(
      PN_MHQ2 == 1 | EN_MHQ2 == 1 | EA_MHQ2 == 1 | PA_MHQ2 == 1 | SA_MHQ2 == 1 ~ 1, # Assign 1 if any of the specified columns are 1
      is.na(ChT_MHQ2) & 
      (PN_MHQ2 == 0 | EN_MHQ2 == 0 | EA_MHQ2 == 0 | PA_MHQ2 == 0 | SA_MHQ2 == 0) ~ 0 # Assign 0 if ChT_MHQ2 is NA and any specified columns are 0
    )
  )
  
#### Merge childhood trauma across MHQ1 and MHQ2 ####

trauma_data$ChT<-NA

trauma_data <- trauma_data %>%
  mutate(
    ChT = case_when(
      ChT_MHQ1 == 1 | ChT_MHQ2 == 1 ~ 1,
      is.na(ChT) & 
      (ChT_MHQ1 == 0 | ChT_MHQ2 == 0 ) ~ 0
      )
    )

save(trauma_data,file="trauma_data_all.Rdata")

#### extract variables for final datasets #####

trauma_final<-trauma_data[,c(1,71:83)]
save(trauma_final,file="Phenotype/CT/trauma_final.Rdata")


  