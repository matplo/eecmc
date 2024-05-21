import uproot
import pandas as pd
from functools import reduce
import re

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
        self.queries = []
        self.groupby_keys = []
        self.n_entries = None
        self.counter = 0
        self.containers = []

    def add(self, file_path, tree_name, prefix, branches=None):
        self.containers.append((file_path, tree_name, prefix, branches))
        self.df = None

    def reset_cuts(self):
        self.min_max_cuts = []
        self.cut_functions = []
        self.queries = []
        self.df = None
        
    def reset_group_by(self):
        self.groupby_keys = []
        self.df = None
        
    def reset(self):
        self.dfs = []
        self.reset_cuts()
        self.reset_group_by()
        self.df = None
        
    def add_group_by(self, *args):
        self.groupby_keys.extend(args)
        self.df = None

    def add_cut(self, prefix, column, lower_bound, upper_bound):
        if column[:len(prefix)+1] != prefix + '_':
            column = prefix + '_' + column
        self.min_max_cuts.append((prefix, column, lower_bound, upper_bound))
        self.df = None

    def add_min_max_cut(self, prefix, column, lower_bound, upper_bound):
        if column[:len(prefix)+1] != prefix + '_':
            column = prefix + '_' + column
        self.min_max_cuts.append((prefix, column, lower_bound, upper_bound))
        self.df = None

    def add_cut_function(self, prefix, func):
        self.cut_functions.append((prefix, func))
        self.df = None

    def add_query(self, prefix, query):
        self.queries.append((prefix, query))

    def rename_columns(self):
        for i, df in enumerate(self.dfs):
            df.rename(columns=lambda x: self.df_prefix[i] + '_' + x if x not in self.groupby_keys else x, inplace=True)

    # def execute(self):
    #     if self.df is None:
    #         self.rename_columns()
    #         self.df = reduce(lambda left,right: pd.merge(left,right,on=self.groupby_keys, how='outer'), self.dfs)
    #         self.df.drop_duplicates(subset=self.groupby_keys, keep='first', inplace=True)
    #         # self.df.dropna(inplace=True)
    #     for column, lower_bound, upper_bound in self.min_max_cuts:
    #         self.df = self.df[(self.df[column] >= lower_bound) & (self.df[column] < upper_bound)]
    #     for func in self.cut_functions:
    #         self.df = func(self.df)
    #     self.df = self.df.groupby(self.groupby_keys)
    #     self.iterator = self.df.__iter__()
    #     self.executed = True

    def load_dfs(self):
        for file_path, tree_name, prefix, branches in self.containers:
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

    def execute(self):
        if self.df is None:
            self.load_dfs()
            self.rename_columns()
            for i, _df in enumerate(self.dfs):
                for _c, _col, _lower_bound, _upper_bound in self.min_max_cuts:
                    if self.df_prefix[i] == _c or _c == '*':
                        _df = _df[(_df[_col] >= _lower_bound) & (_df[_col] < _upper_bound)]
                for _c, _func in self.cut_functions:
                    if self.df_prefix[i] == _c or _c == '*':
                        _df = _func(_df)
                for _c, _query in self.queries:
                    if self.df_prefix[i] == _c or _c == '*':
                        _df = _df.query(_query)
                self.dfs[i] = _df
            self.df = reduce(lambda left,right: pd.merge(left,right,on=self.groupby_keys, how='inner'), self.dfs)
            # self.df.drop_duplicates(subset=self.groupby_keys, keep='first', inplace=True)
            # self.df.drop_duplicates(subset=self.df.columns, keep='first', inplace=True)
            # self.df.drop_duplicates(subset=self.keys, keep='first', inplace=True)
            self.df = self.df.groupby(self.groupby_keys)
        self.iterator = self.df.__iter__()
        self.executed = True

    def _rebuild(self, df):
        _d = {}
        _d['index'] = self.counter
        for k in df.columns:
            _d[k] = df[k].values
        _df = pd.DataFrame(_d)
        # print('rebuild', self.counter)
        # print(_df)
        self.counter += 1
        return _df

    def execute_test2(self):
        if self.df is None:
            self.rename_columns()
            for i, _df in enumerate(self.dfs):
                for _c, _col, _lower_bound, _upper_bound in self.min_max_cuts:
                    if self.df_prefix[i] == _c or _c == '*':
                        _df = _df[(_df[_col] >= _lower_bound) & (_df[_col] < _upper_bound)]
                for _c, _func in self.cut_functions:
                    if self.df_prefix[i] == _c or _c == '*':
                        _df = _func(_df)
                if i > 0:
                    if i == 1:
                        self.df = pd.merge(self.dfs[0], _df, on=self.groupby_keys, how='outer')
                        # print('first merge')
                        # print (self.df)
                    else:
                        self.df = pd.merge(self.df, _df, on=self.groupby_keys, how='outer')
                        # print(f'{i} merge')
                        # print (self.df)
            self._gby = self.df.groupby(self.groupby_keys)
            self.df = self._gby.apply(self._rebuild)
            # print('after groupby')
            # print(self.df)
        self.iterator = self.df.__iter__()
        self.executed = True
    
    def get_unique(self, columns):
        if self.df is None:
            self.execute()
        cols = self.groupby_keys.copy()
        for x in columns:
            if x not in self.df.obj.columns:
                if '*' in x:
                    _ = [cols.append(c) for c in self.df.obj.columns if re.match(x.replace('*', '.*'), c)]
                    continue
                continue
            cols.append(x)
        # keep only unique values in cols
        cols = list(set(cols))
        _df = self.df.apply(lambda x: x[cols].drop_duplicates(cols))
        _df.reset_index(drop=True, inplace=True)
        # return _df[columns]
        return _df
    
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