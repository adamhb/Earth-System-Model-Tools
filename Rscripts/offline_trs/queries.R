#queries
luqfull_raw %>% filter(Status == "broken below", Census == 3) %>% select("TreeID", "Codes", "Status") %>% arrange(TreeID) %>% as.data.frame()
