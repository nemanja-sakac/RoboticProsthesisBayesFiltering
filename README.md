# Bayes Filtering for Robotic Prosthesis EMG Signal Filtering
Nonlinear Bayesian filtering developed for the estimation of an underlying "driving process" (envelope) of a signal. Based on paper and algorithm as demonstrated by Sanger[^1].

The implemented function returns the envelope of the input signal with the possibility of using multiple channels as inputs. Depending on the underlying "driving process", two models of the signal have been implemented:
* Gaussian - Models the signal as amplitude-modulated zero-mean Gaussian noise.
* Laplacian - Models the signal with a Laplace distribution.

The output of the function is an EMG signal envelope which can be used as an input for machine learning (recognition) algorithms. Based on the classification result, the robotic hand prosthesis is able to perform an action initiated by the user of the prosthesis.

![emg9](https://user-images.githubusercontent.com/63115088/153922993-f08e1c69-9346-4d58-bd8d-6ae8e07c9b44.png)


[^1]: Sanger, T.D., "Bayesian Filtering of Myoelectric Signals", Journal of Neurophysiology, vol. 97, pp. 1839-1845, 2007.
