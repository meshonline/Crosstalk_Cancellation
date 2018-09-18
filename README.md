# Crosstalk_Cancellation
C++ class for realtime audio crosstalk cancellation.

This is a C++ version of [https://github.com/Curly-Mo/crosstalk_cancellation](https://github.com/Curly-Mo/crosstalk_cancellation), header only, easy to use.

### How to use
There are only five public methods in the class:<br/>
<br/>
void change_sample_rate(const int _sr);<br/>
void change_speaker_to_speaker(const double _spkr_to_spkr);<br/>
void change_listener_to_speaker(const double _lstnr_to_spkr);<br/>
void change_ear_to_ear(const double _ear_to_ear);<br/>
void process_stereo_channel(std::vector<double>& left, std::vector<double>& right);<br/>
<br/>
If the raw audio sample format is unsigned short, please convert the audio samples to double:<br/>
double new_signal = (raw_signal / 65536.0) - 0.5;<br/>
<br/>
After called process_stereo_channel method, convert the audio samples back to unsigned short:<br/>
unsigned short raw_signal = (new_signal + 0.5) * 65536.0;<br/>
<br/>
Note: process_stereo_channel method always returns the result of the previous input vector, not that of the current input vector, so the length of the output vector may not equals the length of the input vector.

### License
The MIT License (MIT)
