# Crosstalk_Cancellation
C++ class of crosstalk cancellation

This is a port version of [https://github.com/Curly-Mo/crosstalk_cancellation](https://github.com/Curly-Mo/crosstalk_cancellation), I re-implemented it with C++, header only, just include the header to use.

### How to use
There are only five methods in the class:<br/>
void change_sample_rate(const int _sr);<br/>
void change_speaker_to_speaker(const double _spkr_to_spkr);<br/>
void change_listener_to_speaker(const double _lstnr_to_spkr);<br/>
void change_ear_to_ear(const double _ear_to_ear);<br/>
void process_stereo_channel(std::vector<double>& left, std::vector<double>& right);<br/>
<br/>
If the raw audio sample format is unsigned short, please convert it to double, like this:<br/>
double signal = (raw_signal / 65536.0) - 0.5;<br/>
Then call 'process_stereo_channel method', when finished, convert it back to unsigned short:<br/>
unsigned short raw_signal = (signal + 0.5) * 65536.0;

### License
The MIT License (MIT)
