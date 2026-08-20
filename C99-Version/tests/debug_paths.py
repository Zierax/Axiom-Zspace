import os
HERE = os.path.dirname(os.path.abspath(__file__))
print('HERE:', HERE)
tmp = os.path.join(HERE, 'tmp')
print('tmp:', tmp)
lc_path = os.path.join(tmp, 'lc.csv')
print('lc_path:', lc_path)
def to_wsl(p):
    if not p: return ''
    p = p.replace('\\', '/')
    if p.startswith('D:/'): return p.replace('D:/', '/mnt/d/')
    if p.startswith('C:/'): return p.replace('C:/', '/mnt/c/')
    return p
print('to_wsl(lc_path):', to_wsl(lc_path))
