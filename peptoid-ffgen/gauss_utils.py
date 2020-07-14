import os.path
import re
import sys
import argparse
import numpy as np

def check_post_hartree_fock(post_hartree_fock):
    """Check if post-HF method is being used and which one.

    Parameters
    ----------
    post_hartree_fock : str
        String of post-HF level

    Returns
    -------
    using_post_HF : bool
    """
 post_HF_options = {'MP2': 'EUMP2', 'MP3': 'EUMP3'}
    if post_hartree_fock:
        assert post_hartree_fock in post_HF_options.keys(), \
            "{} not available. ".format(post_hartree_fock) +\
            "Must choose from: {}".format(post_HF_options.keys())
    else:
        pass
