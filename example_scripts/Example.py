# Author: Devang Thakkar
# Last Modified by: Devang Thakkar

class Example(Module):
    def __init__(self, module_id, is_docker=False):
        super(Example, self).__init__(module_id, is_docker)

        # Add output keys here if reuired
        self.output_keys        = ["header_lines", "header_chars"]


    def define_input(self):
        self.add_argument("bam",            is_required=True)
        self.add_argument("nr_cpus",        default_value=1)
        self.add_argument("mem",            default_value=5)


    def define_output(self):
        header_lines            = self.generate_unique_file_name("header_lines.txt")
        header_chars            = self.generate_unique_file_name("header_chars.txt")
        self.add_output("header_lines",    header_lines)
        self.add_output("header_chars",    header_chars)


    def define_command(self):
        bam                     = self.get_argument("bam")
        lines                     = self.get_output("header_lines")
        chars                     = self.get_output("header_chars")

        cmd = "bash /script.sh -b {0} -l {1} -c {2}".format(bam, lines, chars)

        # add logging
        cmd += "!LOG3!"

        return cmd

