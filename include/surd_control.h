struct surd_control
{
	static unsigned int DEBUG; //bitmap
	enum parameters
	{
		general = 	0x00000001, // unused
		sanitize =  0x00000002, 
		add = 		0x00000004,
		subtract =	0x00000008,
		times = 	0x00000010,
		equal = 	0x00000020,
		greater = 	0x00000040,
		input = 	0x00000080,
		newton = 	0x00000100,
		// 'all' also supported
	};
};
