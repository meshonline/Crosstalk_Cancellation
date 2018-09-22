//
//  CrossTalkCancellation.hpp
//  CAPlayThrough
//
//  Created by MINGFENWANG on 2018/9/14.
//

#ifndef CrossTalkCancellation_hpp
#define CrossTalkCancellation_hpp

#include <vector>

class CrossTalkCancellation {
public:
    CrossTalkCancellation() {
        sr = 44100;
        spkr_to_spkr = 0.3048;
        lstnr_to_spkr = 0.5588;
        ear_to_ear = 0.215;
        headshadow_filter_coefficients(compute_geometry(), ear_to_ear/2.0);

        // init left queue
        left_queue.resize(4);
        left_queue[0].resize(512, 0.0);
        left_queue[1].resize(512, 0.0);
        left_queue[2].resize(512, 0.0);
        left_queue[3].resize(512, 0.0);

        // init right queue
        right_queue.resize(4);
        right_queue[0].resize(512, 0.0);
        right_queue[1].resize(512, 0.0);
        right_queue[2].resize(512, 0.0);
        right_queue[3].resize(512, 0.0);
    }
    
    void change_sample_rate(const int _sr) {
        sr = _sr;
        headshadow_filter_coefficients(compute_geometry(), ear_to_ear/2.0);
    }
    
    void change_speaker_to_speaker(const double _spkr_to_spkr) {
        spkr_to_spkr = _spkr_to_spkr;
        headshadow_filter_coefficients(compute_geometry(), ear_to_ear/2.0);
    }
    
    void change_listener_to_speaker(const double _lstnr_to_spkr) {
        lstnr_to_spkr = _lstnr_to_spkr;
        headshadow_filter_coefficients(compute_geometry(), ear_to_ear/2.0);
    }
    
    void change_ear_to_ear(const double _ear_to_ear) {
        ear_to_ear = _ear_to_ear;
        headshadow_filter_coefficients(compute_geometry(), ear_to_ear/2.0);
    }
    
    void process_stereo_channel(std::vector<double>& left,
                                std::vector<double>& right) {
        // move old data pieces forward
        left_queue[0] = left_queue[1];
        left_queue[1] = left_queue[2];
        left_queue[2] = left_queue[3];
        // queue new data piece
        left_queue[3] = left;

        // move old data pieces forward
        right_queue[0] = right_queue[1];
        right_queue[1] = right_queue[2];
        right_queue[2] = right_queue[3];
        // queue new data piece
        right_queue[3] = right;

        // concatenate data pieces to one big data piece
        std::vector<double> work_left;
        work_left.insert(work_left.end(), left_queue[0].begin(), left_queue[0].end());
        work_left.insert(work_left.end(), left_queue[1].begin(), left_queue[1].end());
        work_left.insert(work_left.end(), left_queue[2].begin(), left_queue[2].end());
        work_left.insert(work_left.end(), left_queue[3].begin(), left_queue[3].end());

        std::vector<double> work_right;
        work_right.insert(work_right.end(), right_queue[0].begin(), right_queue[0].end());
        work_right.insert(work_right.end(), right_queue[1].begin(), right_queue[1].end());
        work_right.insert(work_right.end(), right_queue[2].begin(), right_queue[2].end());
        work_right.insert(work_right.end(), right_queue[3].begin(), right_queue[3].end());

        // calculate crosstalk cancellation for left channel
        std::vector<double> l_left, l_right;
        l_left.resize(work_left.size(), 0.0);
        l_right.resize(work_left.size(), 0.0);
        cancel_crosstalk(work_left, l_left, l_right);
        
        // calculate crosstalk cancellation for right channel
        std::vector<double> r_right, r_left;
        r_right.resize(work_right.size(), 0.0);
        r_left.resize(work_right.size(), 0.0);
        cancel_crosstalk(work_right, r_right, r_left);
        
        // accumulate crosstalk cancellation to left channel
        add_to_signal(work_left, l_left);
        add_to_signal(work_left, r_left);
        
        // accumulate crosstalk cancellation to right channel
        add_to_signal(work_right, r_right);
        add_to_signal(work_right, l_right);
        
        // store the third data piece as the result
        left.resize(left_queue[2].size());
        memcpy(left.data(),
               work_left.data() + left_queue[0].size() + left_queue[1].size(),
               left_queue[2].size() * sizeof(double));
        right.resize(right_queue[2].size());
        memcpy(right.data(),
               work_right.data() + right_queue[0].size() + right_queue[1].size(),
               right_queue[2].size() * sizeof(double));
    }
private:
    /******************************************************************************
     IIR(Infinite impulse response) filter
     Where:
     b[i] are the feedforward filter coefficients.
     a[i] are the feedback filter coefficients.
     *****************************************************************************/
    struct HeadShadow {
        double b[2];
        double a[2];
    };
    
    
    void cancel_crosstalk(const std::vector<double>& signal,
                          std::vector<double>& ipsilateral,
                          std::vector<double>& contralateral) {
        double c = 343.2;
        double delta_d = abs(d2 - d1);
        double time_delay = delta_d / c;
        double attenuation = d1 / d2;
        // Reference max amplitude
        double ref = max_of_abs(signal);
        recursive_cancel(signal, ipsilateral, contralateral, ref, time_delay, attenuation);
    }
    
    inline void add_to_signal(std::vector<double>& signal,
                       const std::vector<double>& signal2) {
        for (int i=0; i<static_cast<int>(signal.size()); i++) {
            signal[i] += signal2[i];
        }
    }
    
    /******************************************************************************
                ----------------        ----------------
               |   M samples    |      |delta fractional|
     y[n] ---> | integer delay  | ---> |      delay     | ---> y[n - M - delta]
               |                |      |                |
                ----------------        ----------------
     
     The M samples integer delay is implemented as a simple delay line:
     y[n] = x[n - M]
     
     The delta fractional delay is implemented as:
     y[n] = u * x[n] + x[n - 1] - u * y[n - 1] where u = (1 - delta) / (1 + delta)
     
     The formula has initial condition issues, so I use naive method instead.
     *****************************************************************************/
    inline void fractional_delay(std::vector<double>& signal,
                          double time) {
        std::vector<double> temp_signal(signal.size());
        double m = time * sr;
        int m_int = floor(m);
        double m_frac = m - m_int;
        for (int i=0; i<static_cast<int>(temp_signal.size()); i++) {
            int index_low = i - (m_int + 1);
            int index_high = i - m_int;
            double low = (index_low >= 0 ? signal[index_low] : 0.0);
            double high = (index_high >= 0 ? signal[index_high] : 0.0);
            temp_signal[i] = high + (low - high) * m_frac;
        }
        signal = temp_signal;
    }
    
    /******************************************************************************
     A linear filter that achieves zero phase delay by applying an IIR filter to a
     signal twice, once forwards and once backwards. The order of the filter is
     twice the original filter order.
     y[t] = b[0] * x[t] + b[1] * x[t-1] + ... - a[1] * y[t-1] - a[2] * y[t-2] - ...
     To simplify the caculation, I let y[0] = x[0], though this method has initial
     condition issue, I use padding data at both head and tail to avoid the issue.
     *****************************************************************************/
    inline void filtfilt(std::vector<double>& signal) {
        if (signal.size() > 1) {
            std::vector<double> temp_signal;
            // from left to right
            temp_signal = signal;
            for (int i=1; i<static_cast<int>(signal.size()); i++) {
                temp_signal[i] = headshadow.b[0] * signal[i] +
                headshadow.b[1] * signal[i-1] -
                headshadow.a[1] * temp_signal[i-1];
            }
            
            // from right to left
            signal = temp_signal;
            for (int i=static_cast<int>(signal.size())-2; i>=0; i--) {
                signal[i] = headshadow.b[0] * temp_signal[i] +
                headshadow.b[1] * temp_signal[i+1] -
                headshadow.a[1] * signal[i+1];
            }
        }
    }
    
    inline void attenuate(std::vector<double>& signal, const double attenuation) {
        for (int i=0; i<static_cast<int>(signal.size()); i++) {
            signal[i] *= attenuation;
        }
    }
    
    /******************************************************************************
     To speed up calculation and optimize memory, I use loop instead of recursive,
     accumulate crosstalk cancellation to each channel in turn.
     *****************************************************************************/
    void recursive_cancel(const std::vector<double>& signal,
                          std::vector<double>& ipsilateral,
                          std::vector<double>& contralateral,
                          const double ref,
                          const double time,
                          const double attenuation,
                          const double threshold_db = -70.0) {
        std::vector<double> temp_signal = signal;
        bool ping_pong = false;
        double db;
        do {
            // time delay by linear interpolation
            fractional_delay(temp_signal, time);

            // invert the delayed signal
            invert(temp_signal);

            // apply headshadow filter (lowpass based on theta)
            filtfilt(temp_signal);

            // attenuate the low pass filtered delayed signal
            attenuate(temp_signal, attenuation);

            // accumulate signal to either ipsilateral or contralateral
            add_to_signal(ping_pong ? ipsilateral : contralateral,
                          temp_signal);

            // Recurse until rms db is below threshold
            db = 20 * log10(max_of_abs(temp_signal) / ref);

            // flip left and right
            ping_pong = !ping_pong;
        } while (db >= threshold_db);
    }
    
    inline double max_of_abs(const std::vector<double>& signal) {
        double max_abs_x = 0.0;
        for (int i=0; i<static_cast<int>(signal.size()); i++) {
            double abs_x = std::abs(signal[i]);
            if (abs_x > max_abs_x) {
                max_abs_x = abs_x;
            }
        }
        return max_abs_x;
    }
    
    void headshadow_filter_coefficients(const double _theta,
                                        const double r) {
        double theta = _theta + M_PI / 2.0;
        double theta0 = 2.618;
        double alpha_min = 0.5;
        double c = 343.2;
        double w0 = c / r;
        double alpha = 1 + alpha_min / 2.0 + (1.0 - alpha_min / 2.0) * cos(theta * M_PI / theta0);
        
        headshadow.b[0] = (alpha + w0 / sr) / (1 + w0 / sr);
        headshadow.b[1] = (-alpha + w0 / sr) / (1 + w0 / sr);
        headshadow.a[0] = 1;
        headshadow.a[1] = -(1 - w0 / sr) / (1 + w0 / sr);
    }
    
    double compute_geometry() {
        double S = spkr_to_spkr / 2.0;
        double L = lstnr_to_spkr;
        double r = ear_to_ear / 2.0;
        double theta = acos(S / (sqrt(L * L + S * S)));
        double delta_d = r * (M_PI - 2.0 * theta);
        d1 = sqrt(L * L + (S - r) * (S - r));
        d2 = d1 + delta_d;
        
        // angle from center of head to speaker (used for computing headshadow)
        theta = atan(S / L);
        
        return theta;
    }
    
    inline void invert(std::vector<double>& signal) {
        for (int i=0; i<static_cast<int>(signal.size()); i++) {
            signal[i] = -signal[i];
        }
    }
    
    double d1;
    double d2;
    HeadShadow headshadow;
    double spkr_to_spkr;
    double lstnr_to_spkr;
    double ear_to_ear;
    int sr;
    std::vector< std::vector<double> > left_queue;
    std::vector< std::vector<double> > right_queue;
};

#endif /* CrossTalkCancellation_hpp */
