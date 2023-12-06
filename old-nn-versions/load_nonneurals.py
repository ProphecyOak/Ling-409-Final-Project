import json, subprocess
import importlib

def load_pyscript(version_id, version_num):
    base_url = 'https://github.com/ProphecyOak/Ling-409-Final-Project/blob/VERSION_ID/nonneural.py'
    url = base_url.replace('VERSION_ID', version_id)
    subprocess.run(['wget', '-O', 'nonneural.py.json', url])

    with open('nonneural.py.json') as jsonfile:
        nneural_json = json.load(jsonfile)

    with open('nonneural_v{}.py'.format(version_num), 'w') as outfile:
        nneural_raw = nneural_json['payload']['blob']['rawLines']
        outfile.write('\n'.join(nneural_raw))
    
    return nneural_raw

def get_func_names(pymod, raw):
    imported_libs = []

    for line in raw:
        if 'import' in line:
            imported_libs += [mod.replace(' ', '') for mod in line.split('import')[-1].split(',')]
    
    return [name for name in dir(pymod) if name[:2] != '__' and name not in imported_libs]

def load_nneural_map():
    version_ids = [('orig','dc9b2c4ddc94f3064458c12571a80833e1971b78'),
                   ('test_and_dbg','453aeec685f6c04f339b30b43ef70ce11229170e'),
                   ('final','8a9bbb7684530e4eb46f5075eb6006cc35f3d9c6')]
    func_map = {}

    for i,pair in enumerate(version_ids):
        vnum = i+1
        name,sha = pair

        raw_file = load_pyscript(sha, vnum)
        nneural = importlib.import_module('nonneural_v{}'.format(vnum))
        func_names = get_func_names(nneural, raw_file)

        for func_name in func_names:
            if func_name not in func_map.keys():
                func_map[func_name] = {}
            func_map[func_name][name] = getattr(nneural, func_name)
    
    return func_map

