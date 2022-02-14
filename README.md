# Bayes Filtering for Robotic Prosthesis EMG Signal Filtering
<p>Nonlinear Bayesian filtering developed for the estimation of an underlying "driving process" of a signal. Based on paper and algorithm demonstrated in the paper Sanger, T.D., "Bayesian Filtering of Myoelectric Signals", Journal of Neurophysiology, vol. 97, pp. 1839-1845, 2007.<\p>
<p>The function returns the envelope of the input signal. Depending on the underlying "driving process" two models of the signal have been implemented:
* Gaussian - Models the signal as amplitude-modulated zero-mean Gaussian noise.
* Laplacian - Models the signal with a Laplace distribution.
<\p>
<p> The output of the function is an EMG signal envelope that can be used as an input for machine learning (recognition) algorithms. Based on the envelope, the robotic hand prosthesis is able to perform an action wished for by the user using the prosthesis<\p>  
