import uproot
import pandas as pd


class UprootTree:
    def __init__(self):
        self.dfs = []
        self.df = None
        self.grouped = None
        self.group_iterator = None

    def add(self, file_path, tree_name, branches=None):
        file = uproot.open(file_path)
        tree = file[tree_name]
        if branches is None:
            branches = tree.keys()
        df = tree.arrays(branches, library="pd")
        self.dfs.append(df)

    def groupby(self, keys):
        self.df = pd.concat(self.dfs)
        self.grouped = self.df.groupby(keys)
        self.group_iterator = self.grouped.__iter__()

    def next(self):
        return next(self.group_iterator)


# use example

# from uproot_tree import UprootTree
# t = UprootTree()
# t.add("inclusive_10_15_jfch.root", "tn_events_jfch")
# t.add("charm_15_30_jfch.root", "tn_correl_jfch")
# t.groupby("nev")

# use
# t.next()
# in a loop

# try:
#     while True:
#         group_name, group_df = tree.next()
#         # do something with group_name and group_df
# except StopIteration:
#     pass
