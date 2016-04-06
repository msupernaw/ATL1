#include "../../../../Optimization/Optimization2/Optimization.hpp"
#include "../../../../Utilities/IO/StreamedDataFile.hpp"




template<typename T>
class ThetaLog : public atl::ObjectiveFunction<T> {
    std::vector<T> Y;
    std::vector<atl::Variable<T> > X;
    atl::Variable<T> logr0; // = -2.6032947;
    atl::Variable<T> logtheta; // = 0.7625692;
    atl::Variable<T> logK; // = 6.7250075;
    atl::Variable<T> logQ; // = -4.7496015;
    atl::Variable<T> logR; // = -3.1889239;
    
    
    
public:
    
    ThetaLog() {
        
    }
    
    void Initialize() {
        
        
        atl::StreamedDataFile data;
        data.open("ThetaLog.dat");
        
        int size = 0;
        data >> size;

        this->Y.resize(size);
        this->X.resize(size);
        
        
        for (int i = 0; i < X.size(); i++) {
            data >> Y[i];
        }
        
        logr0.SetName("logr0");
        logtheta.SetName("logtheta");
        logK.SetName("logK");
        logQ.SetName("logQ");
        logR.SetName("logR");
        
        /**
         * Register parameters.
         */
        this->RegisterHyperParameter(logr0);
        this->RegisterHyperParameter(logtheta);
        logK.SetBounds(4.6, 7.6);
        this->RegisterHyperParameter(logK);
        this->RegisterHyperParameter(logQ);
        this->RegisterHyperParameter(logR);
        
        /**
         * Register random variables.
         */
        for (int i = 0; i < X.size(); i++) {
            this->RegisterRandomVariable(X[i]);
        }
    }
    
    void step(atl::Variable<T>& jnll, const atl::Variable<T>& x1, const atl::Variable<T>& x2, const atl::Variable<T>& logr0, const atl::Variable<T>& logK, const atl::Variable<T>& logtheta, const atl::Variable<T>& logQ) {
        atl::Variable<T> var = atl::exp(logQ);
        atl::Variable<T> m = (x1 + atl::exp(logr0) * (1.0 - atl::pow(atl::exp(x1) / atl::exp(logK), atl::exp(logtheta))));
        jnll += 0.5 * (atl::log(2.0 * M_PI * var) + atl::pow((x2 - m), 2.0) / var);
    }
    
    void obs(atl::Variable<T>& jnll, const atl::Variable<T>& x, const atl::Variable<T>& logR, int i) {
        atl::Variable<T> var = atl::exp(logR);
        jnll += static_cast<T> (0.5) * (atl::log(static_cast<T> (2.0 * M_PI) * var)+(atl::pow(x - Y[i], static_cast<T> (2.0))) / var);
    }
    
    virtual const atl::Variable<T> Evaluate() {
        
        atl::Variable<T> f = 0.0;
        
        
        int timeSteps = Y.size();
        
        for (int i = 1; i < timeSteps; i++) {
            step(f, X[i - 1], X[i], logr0, logK, logtheta, logQ);
        }
        
        for (int i = 0; i < timeSteps; i++) {
            obs(f, X[i], logR, i);
        }
        return f;
    }
    
    void Report() {
        std::ofstream out;
        out.open("ThetaLog.par");
        out << "logr0 = " << logr0 << "\n";
        out << "logtheta = " << logtheta << "\n";
        out << "logK = " << logK << "\n";
        out << "logQ = " << logQ << "\n";
        out << "logR = " << logR << "\n";
        out << "\nRandom Effects:\n";
        out << "X:\n";
        for (int i = 0; i < X.size(); i++) {
            out << X[i] << " ";
        }
    }
    
};

int main(int argc, char** argv) {
    
    //create the objective function 
    ThetaLog< double> objective_function;
    
    //initialize the objective function
    objective_function.Initialize();
    
    //create an instance of a L-BFGS minimizer
    atl::LBFGS< double> fm;
    
    //set the objective function
    fm.SetObjectiveFunction(&objective_function);
    
    //run the minimizer
    fm.Run();
    
    //dump the results
    objective_function.Report();
    return 0;
}

