#!/usr/bin/env python3
import argparse
import os
# read in inputs and overrride certain tags

notinclude = ['ion-species','ionic-minimize',
                  #'latt-scale','latt-move-scale','coulomb-interaction','coords-type',
                  'ion','climbing','pH','ph',  'autodos',
                  'logfile','pseudos','nimages','max_steps','max-steps','fmax','optimizer','opt',
                  'restart','parallel','safe-mode','hessian', 'step', 'Step',
                  'opt-alpha', 'md-steps', 'econv', 'pdos', 'pDOS', 'lattice-type', 'np',
                  'use_jdftx_ionic', 'jdftx_steps', #Tags below here were added by Cooper
                  'bond-fix', 'bader']

def read_inputs(folder, file = 'inputs'):
        # returns a dictionary of inputs
        with open(os.path.join(folder, file), 'r') as f:
            input_text = f.read()
        tags = {}
        for line in input_text.split('\n'):
            if line == '':
                continue
            if line[0] == '#':
                continue
            words = line.split(' ')
            if words[0] in tags:
                if type(tags[words[0]]) != list:
                    tags[words[0]] = [tags[words[0]]]
                tags[words[0]].append(' '.join(words[1:]))
            else:
                tags[words[0]] = ' '.join(words[1:])
        return tags

def write_inputs(inputs_dic, folder):
    # expects a dictionary of inputs with either strings or lists of strings as vals
    text = ''
    for k,v in inputs_dic.items():
        if type(v) == list:
            for vi in v:
                text += k + ' ' + vi + '\n'
        else:
            text += k + ' ' + v + '\n'
    with open(os.path.join(folder, 'inputs'), 'w') as f:
        f.write(text)

def read_line(line):
        line = line.strip()
        if len(line) == 0 or line[0] == '#': 
            return 0,0,'None'
        linelist = line.split()
        
        defaults = ['default','inputs']
        inputs = read_inputs('./')
        
        if linelist[0] in notinclude:
            # command, value, type
            cmd, val, typ = linelist[0], ' '.join(linelist[1:]), 'script'
            if val in defaults:
                val = inputs[cmd]
            return cmd, val, typ
        else:
            if len(linelist) > 1:
                # command, value, type
                cmd, val, typ = linelist[0], ' '.join(linelist[1:]), 'cmd'
                if val in defaults:
                    val = inputs[cmd]
                return cmd, val, typ
            else:
                return line, '', 'cmd'
        
        return 0,0,'None'

def read_commands(inputs, notinclude = notinclude):
    cmds = []
    script_cmds = {}
#        with open(command_file,'r') as f:
    for line in inputs.split('\n'):
#            conv_logger('Line = ', line)
        cmd, val, typ = read_line(line)
        if typ == 'None':
            continue
        elif typ == 'script':
            if cmd in ['pdos','pDOS']:
                if 'pdos' not in script_cmds:
                    script_cmds['pdos'] = []
                script_cmds['pdos'].append(val)
                continue
            script_cmds[cmd] = val
        elif typ == 'cmd':
            tpl = (cmd, val)
            cmds.append(tpl)

#        cmds += [('core-overlap-check', 'none')]
    cmds = add_dos(cmds, script_cmds)
    return cmds, script_cmds

def main(args):
    with open("inputs", "r") as infile:
        inputs = infile.read()
    inputs_dict = read_inputs('./')
    inputs_dict[args.input[0]] = ' '.join(args.input[1:]) # set the value with the associated key to the text that follows the key
    write_inputs(inputs_dict, './')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='change input file items')
    parser.add_argument('-i', '--input', type=str, nargs='+', help=('list a key and values to change in the input file separated by a space'
                                                         'e.g. -i restart FALSE'))
    args = parser.parse_args()
    main(args)