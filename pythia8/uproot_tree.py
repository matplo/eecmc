import uproot
import pandas as pd
from functools import reduce

# example cut function
# def cut_func(df):
#     return df[(df['column1'] > 0) & (df['column2'] < 1) & (df['column3'] != df['column4'])]
# t = UprootTree()
# t.add("file.root", "tree")
# t.add_cut(cut_func)
# t.execute()


class UprootTree:
    def __init__(self):
        self.dfs = []
        self.df = None
        self.grouped = None
        self.group_iterator = None
        self.min_max_cuts = []
        self.cut_functions = []
        self.groupby_keys = []
        self.df_merged = None
        self.df_copy = None
        self.grouped = None
        self.n_entries = None

    def add(self, file_path, tree_name, branches=None):
        file = uproot.open(file_path)
        tree = file[tree_name]
        if branches is None:
            branches = tree.keys()
        if self.n_entries is None:
            df = tree.arrays(branches, library="pd")
        else:
            df = tree.arrays(branches, library="pd", entry_stop=self.n_entries)
        self.dfs.append(df)
        self.df_merged = None

    def reset_cuts(self):
        self.min_max_cuts = []
        self.cut_functions = []
        self.grouped = None
        
    def reset_group_by(self):
        self.groupby_keys = []
        self.grouped = None
        
    def reset(self):
        self.reset_cuts()
        self.reset_group_by()
        self.df_merged = True
        
    def add_group_by(self, *args):
        self.groupby_keys.extend(args)
        self.grouped = None

    def add_cut(self, column, lower_bound, upper_bound):
        self.min_max_cuts.append((column, lower_bound, upper_bound))
        self.grouped = None

    def add_min_max_cut(self, column, lower_bound, upper_bound):
        self.min_max_cuts.append((column, lower_bound, upper_bound))
        self.grouped = None

    def add_cut_function(self, func):
        self.cut_functions.append(func)
        self.grouped = None

    def execute(self):
        if self.df_merged is None:
            # self.df_merged = pd.concat(self.dfs, join='outer')
            self.df_merged = reduce(lambda left,right: pd.merge(left,right,on=self.groupby_keys, how='outer'), self.dfs)
            # self.df_merged.fillna(0, inplace=True) # fill missing values with zero's
            self.df_merged.dropna(inplace=True) # drop rows with missing values
            # print(self.df_merged)            
        for column, lower_bound, upper_bound in self.min_max_cuts:
            self.df_copy = self.df_merged[(self.df_merged[column] >= lower_bound) & (self.df_merged[column] < upper_bound)]
        for func in self.cut_functions:
            self.df_copy = func(self.df_copy)
        self.df_copy.reset_index(drop=True, inplace=True)
        self.grouped = self.df_copy.groupby(self.groupby_keys)
        # self.df = self.grouped.reset_index(drop=True)
        self.df = self.grouped.apply(lambda x: x).reset_index(drop=True)
        self.group_iterator = self.grouped.__iter__()
        self.executed = True

    def reset_iterator(self):
        self.group_iterator = self.grouped.__iter__()
        
    def next(self):
        if self.grouped is None:
            self.df_copy = None
            self.execute()
        return next(self.group_iterator)

    def get_df(self):
        if self.grouped is None:
            self.df_copy = None
            self.execute()
        return self.df