class wFunction():
    def __init__(self, config, section, function):

        self.config = config.copy()
        self.section = section
        self.function = function

        self._read_config()

        pass

    # Override to prepare the function
    # Calculate and to add the remaining parameters for the function
    def preparation(self):
        pass

    # Override to adjust the output of the function
    # Modify the output if needed
    def output(self, output):
        pass

    # Do not override
    def _read_config(self):
        import inspect as ins

        par_list = dict()

        # Inspect the function and get all the parameter required
        arg = ins.signature(self.function).parameters
        for key in arg.keys():
            if self.config.has_option(self.section, key):
                par_list[key] = self.config[self.section][key]

        self._par_list = par_list

    # Do not override
    def run(self):
        self.preparation()
        output = self.function(**self._par_list)
        self.output(output)
