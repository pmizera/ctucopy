/*
 CtuCopy 3 - Universal feature extractor and speech enhancer.
 Copyright 2012 Petr Fousek, FEE CTU Prague

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/

/// mu-law / A-law to linear PCM convertor

void alaw2lin (char &a, short int &out, int mode = 1) {
	// mode: 0=mu-law, 1=A-law
	
	int chord;       // 3 bit pcm chord
	int step;        // 4 bit pcm step
	int mag;
	
	int sgn = (~(a >> 7)) & 1; // sign bit
	
	if (mode == 0)
		{
			chord = (~(a >> 4)) & 7;
			step = (~(a)) & 0xf;
			mag = (((2 * step) + 33) << chord) - 33;
		}
	else
		{ 
			// ALAW linear number is 2x CCITT value in order to match mu-Law
			chord = ((a ^ 0x55) >> 4) & 7;
			step  = ((a ^ 0x55)) & 0xf;
			mag   = (step << 1) + 1;
			if (chord > 0)
				mag += 32;
			else
				chord = 1;
			mag = mag << chord;
		}
	
		out = ((1 - (2 * sgn)) * mag)  & 0xffff;

		out = (out << 2) & 0xffff; // 2x amplification

		if ((out & 0x8000) != 0 ) out = out - 65536; 
}















