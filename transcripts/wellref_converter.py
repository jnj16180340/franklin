#!/usr/bin/python3
import itertools
import string
import math
import argparse

def num2let(arg):
    '''Converts integer [0...25] to character [A...Z]. Fails more noticeably than chr()!'''
    return dict(zip(range(26), string.ascii_uppercase))[arg]

def let2num(arg):
    '''Converts character [A...Z] to integer [0...25]. Fails more noticeably than ord()!'''
    return dict(zip(string.ascii_uppercase, range(26)))[arg]

def human2robot_column(human_column_string):
    '''Converts integer: Robot numbering (0 based) -> human numbering (1 based)'''
    return(int(human_column_string)-1)

def robot2human_column(robot_column_integer):
    '''Converts integer: Human numbering (1 based) -> robot numbering (0 based)'''
    return(str(robot_column_integer+1))

# TODO: For tests
#assert human2robot_column(robot2human_column(0)) ==  0, 'bad'
#assert robot2human_column(human2robot_column('0')) == '0', 'bad'

def robot2human_row(robot_row_int,base=26):
    '''Converts integer: Base 10 -> character string representation of base N (default 26)
    See https://en.wikipedia.org/wiki/Positional_notation#Base_conversion
    '''
    (q,r) = divmod(robot_row_int,base)
    if q:
        return robot2human_row(q-1) + num2let(r)
    else:
        return num2let(r)

def human2robot_row(human_row_string):
    '''converts string rep base 26/1 start -> digit[0...9] rep base 10/0 start'''
    result = -1
    exponent = 0
    for k in reversed(human_row_string): 
        result += (let2num(k)+1) * 26**exponent
        exponent+=1
    return(result)

#TODO: for tests
#for k in [human2robot_row(robot2human_row(k)) == k for k in range(500)]:
#    assert k, "bad"
    
for k in [human2robot_column(robot2human_column(k)) == k for k in range(500)]:
    assert k, "bad"

def robotize_to_tuple(human_well_string):
    ''' Converts string representation of a single well into tuple(int,int) representation of well
    e.g. "AB13" -> (28,13)
    '''
    filter(str.isdigit,human_well_string)
    row = ''.join(filter(str.isalpha,human_well_string))
    column = int(''.join(filter(str.isdigit,human_well_string)))
    return human2robot_row(row),human2robot_column(column)

def humanize_from_tuple(robot_well_tuple):
    '''Converts tuple(int,int) rep of a single well into string representation of well
    e.g. (28,13) -> "AB13"
    '''
    return robot2human_row(robot_well_tuple[0])+robot2human_column(robot_well_tuple[1])

def robot_tuple2int(robot_tuple_well, platedims):
    '''Converts tuple(int,int) representation of well into column-based int representation of well
    e.g. 3 -> (1,1) for standard 4-well plate
    '''
    platecolumns,platecols = platedims
    r,c = robot_tuple_well
    return r*platecols+c

def robot_int2tuple(robot_int_well, platedims):
    '''Converts column-based int representation of well into tuple(int,int) representation of well
    e.g. (1,1) -> 3 for standard 4-well plate
    '''
    platecolumns,platecols = platedims
    return divmod(robot_int_well,platecols)

# TODO: for tests
#assert humanize_from_tuple(robotize_to_tuple('A1')) == 'A1', 'bad'
#assert robotize_to_tuple(humanize_from_tuple((11,3))) == (11,3), 'bad'
#assert robot_int2tuple(robot_tuple2int((1,0),(8,12)),(8,12)) == (1,0), 'bad'
#assert robot_tuple2int(robot_int2tuple(99,(8,12)),(8,12)) == 99, 'bad'

def robotize(well, platedims):
    '''Converts human string representation of well into robot int representation of well
    e.g. "A1" -> 0
    '''
    return robot_tuple2int(robotize_to_tuple(well),platedims)

def humanize(well, platedims):
    '''Converts robot int representation of well into human string representation of well
    e.g. "A1" -> 0
    '''
    return humanize_from_tuple(robot_int2tuple(well,platedims))

# TODO: for tests
#assert humanize(robotize('A2',(8,12)),(8,12)) == 'A2', 'bad'
#assert robotize(humanize(99,(8,12)),(8,12)) == 99, 'bad'


def main():
    '''Wraps robotize() and humanize() in commandline interface'''
    
    # Common well-plate sizes
    common_plate_sizes = (6, 24, 96, 384, 1536, 3456, 9600) # according to wikipedia!
    common_plate_proportion = (2,3)
    generate_common_plate_dimensions = lambda size: tuple([k * int(math.sqrt(size/6)) for k in common_plate_proportion])
    common_plate_dimensions = dict(zip((''.join(str(k)+'well') for k in common_plate_sizes), map(generate_common_plate_dimensions, common_plate_sizes)))
    
    parser = argparse.ArgumentParser(description='well ref human <-> robot conversion')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--robotize', help='translate well refs human -> robot', dest='do_robotize', default=False, action='store_true')
    group.add_argument('--humanize', help='translate well refs robot -> human', dest='do_humanize', default=False, action='store_true')
    parser.add_argument('wells', nargs='+', help='well references, e.g.: (human) A1 B2 J6 AA70 ABZ900... LetterNumber (robot) 10 5 1000... Number')
    parser.add_argument('--checkbounds', action="store_true", dest="checkbounds", default=False,help='check if wells are outside the dimensions of the plate.')
    parser.add_argument('--platetype', action="store", dest="platetype", type=str,required=False,default='96well', help='plate type: name of common type, e.g. 96well, or dimensions RxC e.g. 8x12')
    args = parser.parse_args()

    if args.platetype in common_plate_dimensions:
        plate_dimensions = common_plate_dimensions[args.platetype]
    else:
        # TODO: plate dimensions could be parsed/validated more robustly
        plate_dimensions = tuple(int(k) for k in args.platetype.split('x'))
    
    operation = 'robotize' if args.do_robotize else 'humanize' if args.do_humanize else None
    
    # TODO: well refs could be validated much more robustly
    assert all(map(str.isalnum,args.wells)), "well references should contain only A...Z and 0...9"
    
    if args.do_robotize:
        # quick and dirty check to see if each well contains both letters and numbers 
        # TODO make more robust: Check for LetterNumber format, use regex?
        assert all(map(lambda k: all((k.isalnum(),not(any((k.isnumeric(),k.isalpha()))))), args.wells)), "human-formatted wells should contain both letters and integers"
        output = (str(robotize(well.upper(), plate_dimensions)) for well in args.wells)
        # TODO this wouldn't work if args.wells were an iterator
        if args.checkbounds:
            assert max((robotize(well.upper(), plate_dimensions) for well in args.wells)) < plate_dimensions[0]*plate_dimensions[1], 'well is outside of plate'
        print(' '.join(output))
        
    if args.do_humanize:
        assert all(map(str.isnumeric,args.wells)), "well references should contain only integers"
        # TODO this wouldn't work if args.wells were an iterator
        if args.checkbounds:
            assert max((int(k) for k in args.wells)) < plate_dimensions[0]*plate_dimensions[1], 'well is outside of plate'
        output = (humanize(int(well), plate_dimensions) for well in args.wells)
        print(' '.join(output))
    
if __name__ == '__main__': 
    main()