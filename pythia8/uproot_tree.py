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
        self.df_prefix = []
        self.df = None
        self.iterator = None
        self.min_max_cuts = []
        self.cut_functions = []
        self.groupby_keys = []
        self.n_entries = None

    def add(self, file_path, tree_name, prefix, branches=None):
        file = uproot.open(file_path)
        tree = file[tree_name]
        if branches is None:
            branches = tree.keys()
        if self.n_entries is None:
            df = tree.arrays(branches, library="pd")
        else:
            df = tree.arrays(branches, library="pd", entry_stop=self.n_entries)
        self.dfs.append(df)
        self.df_prefix.append(prefix)
        self.df = None

    def reset_cuts(self):
        self.min_max_cuts = []
        self.cut_functions = []
        self.df = None
        
    def reset_group_by(self):
        self.groupby_keys = []
        self.df = None
        
    def reset(self):
        self.reset_cuts()
        self.reset_group_by()
        self.df = True
        
    def add_group_by(self, *args):
        self.groupby_keys.extend(args)
        self.df = None

    def add_cut(self, column, lower_bound, upper_bound):
        self.min_max_cuts.append((column, lower_bound, upper_bound))
        self.df = None

    def add_min_max_cut(self, column, lower_bound, upper_bound):
        self.min_max_cuts.append((column, lower_bound, upper_bound))
        self.df = None

    def add_cut_function(self, func):
        self.cut_functions.append(func)
        self.df = None

    def rename_columns(self):
        for i, df in enumerate(self.dfs):
            df.rename(columns=lambda x: self.df_prefix[i] + '_' + x if x not in self.groupby_keys else x, inplace=True)

    def execute(self):
        if self.df is None:
            self.rename_columns()
            self.df = reduce(lambda left,right: pd.merge(left,right,on=self.groupby_keys, how='outer'), self.dfs)
            # self.df.dropna(inplace=True)
        for column, lower_bound, upper_bound in self.min_max_cuts:
            self.df = self.df[(self.df[column] >= lower_bound) & (self.df[column] < upper_bound)]
        for func in self.cut_functions:
            self.df = func(self.df)
        # self.df.reset_index(drop=True, inplace=True)
        #self.df = self.df.groupby(self.groupby_keys).apply(lambda x: x).reset_index(drop=True)
        self.df = self.df.groupby(self.groupby_keys)
        self.iterator = self.df.__iter__()
        self.executed = True

    def reset_iterator(self):
        self.iterator = self.df.__iter__()
        
    def next(self):
        if self.df is None:
            self.execute()
        return next(self.iterator)

    def get_df(self):
        if self.df is None:
            self.execute()
        return self.df