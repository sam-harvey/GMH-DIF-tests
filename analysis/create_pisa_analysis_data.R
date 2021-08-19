#Creates the dataset used in the section 5.2 PISA analysis
#This consists of all Australian students who attempted booklet 1
source('libraries.R')

df_scored_cogs = read_csv('data/raw-responses/scored-cognitive-items.csv')
df_student = read_csv('data/responses/student-questionnaire-responses.csv')
df_booklet_questions = read_csv('data/booklet/booklet-items.csv')
df_student_data_dict = read_csv('data/data-dictionaries/student-questionnaire-data-dictionary.csv')

#Join the student's gender and code it
df_aus_students = df_scored_cogs %>% 
  filter(CNT == "AUS",
         BOOKID == 1) %>% 
  inner_join(df_student %>% 
               select(STIDSTD, CNT, ST04Q01)) %>% 
  select(STIDSTD, ST04Q01, starts_with('P')) 

df_aus_students = df_aus_students %>% 
  inner_join(df_student_data_dict %>% 
               filter(VARIABLE_NAME=="ST04Q01") %>% 
               mutate(VALUE = as.numeric(VALUE)) %>% 
               select(VALUE, Gender = VALUE_LABEL), by = c('ST04Q01' = 'VALUE')) %>% 
  select(STIDSTD, Gender, everything())

#Get all science questions in the first book
df_booklet_questions = df_booklet_questions %>% 
  filter(BOOKID == 1) %>% 
  filter(str_detect(COGNITIVE_ITEM_CODE, 'PS'))

test_cols = intersect(colnames(df_aus_students), df_booklet_questions$COGNITIVE_ITEM_CODE)

#Select students to exclude
#We exclude any who did not respond to at least one of the science questions in booklet 1
#This excludes 13 of 1145 total students, leaving 1135 remaining
df_aus_students = df_aus_students %>% 
  select(STIDSTD, Gender, test_cols) %>% 
  filter(if_all(starts_with('PS'), ~. != 7))

write_csv(df_aus_students, 'data/pisa_aus_science_book_1_responses.csv')