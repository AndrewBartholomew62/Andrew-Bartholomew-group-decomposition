struct bigreal_control
{
	static unsigned int DEBUG; //bitmap
	enum parameters
	{
		general = 	0x00000001, // unused
		sanitize =  0x00000002, 
		add = 		0x00000004,
		subtract = 	0x00000008,
		multiply = 	0x00000010,
		divide = 	0x00000020,
		remainder = 0x00000020,
		equal =  	0x00000040,
		greater =  	0x00000080,
	    output =  	0x00000100,
		input =  	0x00000200,
		sum =  		0x00000400,
		diff =  	0x00000800,
		num_len =  	0x00001000,
		r_shift =  	0x00002000,
		l_shift =  	0x00004000,
		bool_conv =	0x00008000,
		gcd =  		0x00010000,
		carry =		0x00020000,
		// 'all' also supported
	};
};
