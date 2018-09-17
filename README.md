# Crosstalk_Cancellation
C++ class of crosstalk cancellation

This is a port version of [https://github.com/Curly-Mo/crosstalk_cancellation](https://github.com/Curly-Mo/crosstalk_cancellation), I re-implemented it with C++, header only.

### How to use
There are only five methods in the class:<br/>
void change_sample_rate(const int _sr);<br/>
void change_speaker_to_speaker(const double _spkr_to_spkr);<br/>
void change_listener_to_speaker(const double _lstnr_to_spkr);<br/>
void change_ear_to_ear(const double _ear_to_ear);<br/>
void process_stereo_channel(std::vector<double>& left, std::vector<double>& right);<br/>
<br/>
If the raw audio sample format is unsigned short, please convert the audio data to double:<br/>
double signal = (raw_signal / 65536.0) - 0.5;<br/>
After called process_stereo_channel, convert the audio data back to unsigned short:<br/>
unsigned short raw_signal = (signal + 0.5) * 65536.0;<br/>
<br/>
Note: The length of the output vector may not equals the length of the input vector.

### License
The MIT License (MIT)
