#Programmer: Nicolas Dullius Mallmann

# TO DO
#   Implement the "=" operator to accept non-string type;

import sys
import json as js
import collections as col

#Some packages depend on the python version
if sys.version_info.major == 2:
    import ConfigParser as cp
elif sys.version_info.major == 3:
    import configparser as cp

class Config(cp.SafeConfigParser):
    """
        Class Config inherits ConfigParser (configparser if python 3).

        Used to read and manipulate a configuration file.

        Configuration File Nomenclature:
            Section = Name of the section where the option will be written.
            Option  = Variable identifier.
            Value   = Value held by the option.

        Useful methods of Config:
            add_section = Add a section.
            copy        = Create a copy of the Config objet.
            get         = Get the value of an option inside a section.
            get_all     = Get every section, option, and value inside an OrderedDict.
            items       = Return a list with the options and values of a section. 
            set         = Set an option with a value (inside a section).
            write       = Create a configuration file with all the sections, options, and values.

        Ex:
            > c1 = Config()
            > c1.add_section("Coord")
            > c1.set("coord", "x", 10)
            > c1.set("coord", "y", -7)
            > c1.set("coord", "y", -3)
            > c1.set("coord", "type", "Cartesian")
            > c2 = c1.copy()
            > c2.write("./config.ini")
            > y = c2.get("coord", "y")
    """

    #Initialization
    def __init__(self, config_file_path = None):

        #Use the father initialization
        cp.SafeConfigParser.__init__(self)

        #Read the configuration file
        if config_file_path != None:
            self.read(config_file_path)

    #Creates a copy
    def copy(self):
        """
        Create a copy of the original object.
        """

        #Create an empty configuration
        config_copy = Config()

        #Read all the parameters from the original configuration
        all_param = self.get_all()

        #Set all the parameters into the empty configuration
        for section in all_param.keys():
            config_copy.add_section(section)

            for option in all_param[section].keys():
                value = self.get(section, option)
                config_copy.set(section, option, value)

        #Return the configuration copy
        return config_copy

    #Override the set function
    def set(self, section, option, value=None):

        #Use the father's set function
        cp.SafeConfigParser.set(self, section, option, str(js.dumps(value)))

    #Override the get function
    if sys.version_info.major == 2:
        def get(self, section, option, raw=False, vars=None):
    
            #Get the original value in string format
            value = cp.SafeConfigParser.get(self, section, option, raw = raw, vars = vars)

            #Return the converted value
            return js.loads(value)
    elif sys.version_info.major == 3:
        def get(self, section, option, raw=False, vars=None, fallback = None):
            
            #Get the original value in string format
            value = cp.SafeConfigParser.get(self, section, option, raw = raw, vars = vars, fallback = fallback)

            #Return the converted value
            return js.loads(value)


    #Override the items function
    def items(self, section, raw=False, vars=None):

        #Get the list in its original format
        list_items = cp.SafeConfigParser.items(self, section = section, raw = raw, vars = vars)

        #Create new list
        new_list = []

        #Convert the value to its type and add a tuple to the new list
        for i in range(len(list_items)):
            item = [list_items[i][0], js.loads(list_items[i][1])]
            new_list.append(tuple(item))

        #Return the list with the new format
        return tuple(new_list)


    #Override the write function
    def write(self, file_path = "./default.ini"):

        #Open a file to write the ini values
        file_handler = open(file_path, "w")

        #Use the father class' original write function with the file handler
        cp.SafeConfigParser.write(self, file_handler)

        #Close the file
        file_handler.close()

    #Return a dictionary with all the values
    def get_all(self):
        """
        Returns an ordered dictionary with every section, option and value of the object.
        """

        #Create the main (ordered) dictionary
        main_dict = col.OrderedDict()

        #Save the section dictionaries inside the main dictionary
        for section in self.sections():
            sect_dict = col.OrderedDict()
            for item in self.items(section, raw = False, vars = None):
                sect_dict[item[0]] = item[1]
            main_dict[section] = sect_dict

        #Return the dictionary
        return main_dict

    #Override the special printable function
    def __str__(self):

        #Initial line
        ret = "\nConfiguration File\n------------------\n"

        #Add the rest of the string
        for section in self.sections():
            ret += " %s:\n" % (section)
            for item in cp.SafeConfigParser.items(self, section, raw = False, vars = None):
                ret += "  -> %s: %s\n" % (item[0], item[1])

        #Return the completed string 
        return ret

    #Override the special representation function
    def __repr__(self):

        #Initial line
        ret = "\nConfiguration File\n------------------\n"

        #Add the rest of the string
        for section in self.sections():
            ret += " %s:\n" % (section)
            for item in cp.SafeConfigParser.items(self, section, raw = False, vars = None):
                ret += "  -> %s: %s\n" % (item[0], item[1])

        #Return the completed string 
        return ret
