import os


def make_directory(target_dir):
    '''
    Make target directory if it does not exist.
    '''
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)


def make_parent_directory(file_name):
    '''
    Given a file name, make the parent directory if it does not exist using make_directory.

    For example, given /some/where/foo.bar make the folder /some/where.
    '''
    file_name = os.path.abspath(file_name)

    parent_dir = os.path.dirname(file_name)

    make_directory(parent_dir)
