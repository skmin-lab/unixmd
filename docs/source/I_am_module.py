class Test_class(object):
    """Test class that sums strings

    :param str a: first string
    :param str b: second string
    """

    def __init__(self, a='I am', b='module'):
        

        self.a = a
        self.b = b

    def string_sum(self):
        """Actual method that sum strings
        
        :returns: string result of a+b
        
        """
        c = self.a + self.b

        return c
