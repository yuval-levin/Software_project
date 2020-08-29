

//TODO: is u1 int? long? double?

void computeS(double* u1,int* s,int n)
{
	int i;
	for(i = 0;i < n;i++) s[i] = u1[i] >= 0 ? 1 : -1;
}

