#ifndef ILOCPLEX_WINDOWS_COMPAT_H_
#define ILOCPLEX_WINDOWS_COMPAT_H_

// ============================================================================
// Camada de Compatibilidade IBM ILOG CPLEX Concert Technology para Windows/GCC
// Permite compilar e executar PL.cpp diretamente no Windows/VS Code usando g++
// mantendo 100% da sintaxe oficial do CPLEX (IloEnv, IloModel, IloCplex, etc.).
// ============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <ctime>

#define ILOSTLBEGIN using namespace std;

const double IloInfinity = 1e30;

enum IloNumVarType {
    ILOFLOAT = 0,
    ILOINT   = 1,
    ILOBOOL  = 2
};

class IloException : public std::exception {
    std::string msg_;
public:
    IloException(const std::string& m = "CPLEX Exception") : msg_(m) {}
    virtual const char* what() const noexcept override { return msg_.c_str(); }
    friend std::ostream& operator<<(std::ostream& os, const IloException& e) {
        return os << e.msg_;
    }
};

struct IloVarData {
    double lb;
    double ub;
    IloNumVarType type;
    std::string name;
};

struct IloRangeData {
    double lhs;
    double rhs;
    std::map<int, double> terms;
    std::string name;
};

struct IloEnvImpl {
    std::vector<IloVarData> vars;
    std::ostringstream nullStream;
};

class IloEnv {
public:
    IloEnvImpl* impl_;
    IloEnv() : impl_(new IloEnvImpl()) {}
    std::ostream& getNullStream() { return impl_->nullStream; }
    std::ostream& error() { return std::cerr; }
    void end() {}
};

class IloNumVar {
public:
    IloEnvImpl* env_;
    int id_;

    IloNumVar() : env_(nullptr), id_(-1) {}
    IloNumVar(IloEnv env, double lb = 0.0, double ub = IloInfinity, IloNumVarType type = ILOFLOAT, const char* name = "")
        : env_(env.impl_) {
        id_ = (int)env_->vars.size();
        IloVarData d;
        d.lb = lb;
        d.ub = ub;
        d.type = type;
        d.name = name ? name : ("x_" + std::to_string(id_));
        env_->vars.push_back(d);
    }

    void setName(const char* name) {
        if (env_ && id_ >= 0 && id_ < (int)env_->vars.size() && name) {
            env_->vars[id_].name = name;
        }
    }
};

class IloNumVarArray {
    std::vector<IloNumVar> arr_;
public:
    IloNumVarArray() {}
    IloNumVarArray(IloEnv, int n = 0) : arr_(n) {}
    void add(const IloNumVar& v) { arr_.push_back(v); }
    IloNumVar& operator[](int i) { return arr_[i]; }
    const IloNumVar& operator[](int i) const { return arr_[i]; }
    int getSize() const { return (int)arr_.size(); }
};

template <typename T>
class IloArray {
    std::vector<T> arr_;
public:
    IloArray() {}
    IloArray(IloEnv, int n = 0) : arr_(n) {}
    T& operator[](int i) { return arr_[i]; }
    const T& operator[](int i) const { return arr_[i]; }
    int getSize() const { return (int)arr_.size(); }
};

class IloExpr {
public:
    std::map<int, double> terms_;
    double constant_;

    IloExpr() : constant_(0.0) {}
    IloExpr(IloEnv) : constant_(0.0) {}
    IloExpr(const IloNumVar& v) : constant_(0.0) {
        if (v.id_ >= 0) terms_[v.id_] += 1.0;
    }

    void clear() {
        terms_.clear();
        constant_ = 0.0;
    }
    void end() { clear(); }

    IloExpr& operator+=(const IloExpr& other) {
        for (auto& kv : other.terms_) terms_[kv.first] += kv.second;
        constant_ += other.constant_;
        return *this;
    }
    IloExpr& operator+=(const IloNumVar& v) {
        if (v.id_ >= 0) terms_[v.id_] += 1.0;
        return *this;
    }
    IloExpr& operator-=(const IloExpr& other) {
        for (auto& kv : other.terms_) terms_[kv.first] -= kv.second;
        constant_ -= other.constant_;
        return *this;
    }
    IloExpr& operator-=(const IloNumVar& v) {
        if (v.id_ >= 0) terms_[v.id_] -= 1.0;
        return *this;
    }
};

inline IloExpr operator*(double coef, const IloNumVar& v) {
    IloExpr e;
    if (v.id_ >= 0) e.terms_[v.id_] = coef;
    return e;
}
inline IloExpr operator*(const IloNumVar& v, double coef) {
    return coef * v;
}
inline IloExpr operator*(int coef, const IloNumVar& v) {
    return (double)coef * v;
}
inline IloExpr operator*(const IloNumVar& v, int coef) {
    return (double)coef * v;
}
inline IloExpr operator+(const IloExpr& a, const IloExpr& b) {
    IloExpr r = a;
    r += b;
    return r;
}
inline IloExpr operator+(const IloExpr& a, const IloNumVar& b) {
    IloExpr r = a;
    r += b;
    return r;
}
inline IloExpr operator+(const IloNumVar& a, const IloNumVar& b) {
    IloExpr r(a);
    r += b;
    return r;
}
inline IloExpr operator-(const IloExpr& a, const IloExpr& b) {
    IloExpr r = a;
    r -= b;
    return r;
}
inline IloExpr operator-(const IloExpr& a, const IloNumVar& b) {
    IloExpr r = a;
    r -= b;
    return r;
}
inline IloExpr operator-(const IloNumVar& a, const IloNumVar& b) {
    IloExpr r(a);
    r -= b;
    return r;
}

class IloRange {
public:
    IloRangeData data_;
    IloRange() {}
    IloRange(IloEnv, double lhs, const IloExpr& exp, double rhs, const char* name = "") {
        data_.lhs = lhs - exp.constant_;
        data_.rhs = rhs - exp.constant_;
        data_.terms = exp.terms_;
        data_.name = name ? name : "";
    }
    IloRange(IloEnv, const IloExpr& exp, double rhs, const char* name = "") {
        data_.lhs = -IloInfinity;
        data_.rhs = rhs - exp.constant_;
        data_.terms = exp.terms_;
        data_.name = name ? name : "";
    }
};

class IloRangeArray {
public:
    std::vector<IloRange> ranges_;
    IloRangeArray() {}
    IloRangeArray(IloEnv) {}
    void add(const IloRange& r) { ranges_.push_back(r); }
    void clear() { ranges_.clear(); }
};

class IloObjective {
public:
    bool isMinimize_;
    std::map<int, double> terms_;
    IloObjective(bool min = true) : isMinimize_(min) {}
    IloObjective(IloEnv, const IloExpr& e, bool min = true) : isMinimize_(min), terms_(e.terms_) {}
    void setLinearCoef(const IloNumVar& v, double c) {
        if (v.id_ >= 0) terms_[v.id_] = c;
    }
};

inline IloObjective IloMinimize(IloEnv env, const IloExpr& e = IloExpr()) {
    return IloObjective(env, e, true);
}
inline IloObjective IloMaximize(IloEnv env, const IloExpr& e = IloExpr()) {
    return IloObjective(env, e, false);
}

class IloModel {
public:
    IloEnvImpl* env_;
    IloObjective obj_;
    std::vector<IloRangeData> ranges_;

    IloModel(IloEnv env) : env_(env.impl_), obj_(true) {}
    void add(const IloObjective& o) { obj_ = o; }
    void add(const IloRange& r) { ranges_.push_back(r.data_); }
    void add(const IloRangeArray& ra) {
        for (const auto& r : ra.ranges_) ranges_.push_back(r.data_);
    }
};

class IloAlgorithm {
public:
    enum Status {
        Unknown    = 0,
        Optimal    = 1,
        Feasible   = 2,
        Infeasible = 3
    };
};

class IloCplex {
    IloEnvImpl* env_;
    IloModel* model_;
    double timeLimit_;
    IloAlgorithm::Status status_;
    double objVal_;
    double gap_;
    std::vector<double> sol_;

public:
    struct Param {
        enum IntParam { TimeLimit = 1 };
    };

    IloCplex(IloEnv env)
        : env_(env.impl_), model_(nullptr), timeLimit_(180.0),
          status_(IloAlgorithm::Unknown), objVal_(0.0), gap_(0.0) {}

    void extract(IloModel& m) { model_ = &m; }
    void exportModel(const char*) {}
    void setOut(std::ostream&) {}
    void setWarning(std::ostream&) {}
    void setParam(int, double val) { timeLimit_ = val; }

    double getCplexTime() const {
        return (double)clock() / CLOCKS_PER_SEC;
    }

    IloAlgorithm::Status getStatus() const { return status_; }
    double getObjValue() const { return objVal_; }
    double getMIPRelativeGap() const { return gap_; }
    double getValue(const IloNumVar& v) const {
        if (v.id_ >= 0 && v.id_ < (int)sol_.size()) return sol_[v.id_];
        return 0.0;
    }

    bool solve() {
        if (!model_ || !env_) return false;
        int nVars = (int)env_->vars.size();
        int nCons = (int)model_->ranges_.size();

        const char* inFile = "_cplex_bridge_in.txt";
        const char* outFile = "_cplex_bridge_out.txt";
        const char* pyFile = "_cplex_bridge_solver.py";

        {
            std::ofstream f(inFile);
            f.precision(15);
            f << nVars << " " << nCons << " " << (model_->obj_.isMinimize_ ? 1 : -1) << " " << timeLimit_ << "\n";
            for (int i = 0; i < nVars; i++) {
                double c = 0.0;
                auto it = model_->obj_.terms_.find(i);
                if (it != model_->obj_.terms_.end()) c = it->second;
                f << env_->vars[i].lb << " " << env_->vars[i].ub << " " << (int)env_->vars[i].type << " " << c << "\n";
            }
            for (int r = 0; r < nCons; r++) {
                const auto& rg = model_->ranges_[r];
                f << rg.lhs << " " << rg.rhs << " " << rg.terms.size();
                for (const auto& kv : rg.terms) {
                    f << " " << kv.first << " " << kv.second;
                }
                f << "\n";
            }
        }

        {
            std::ofstream py(pyFile);
            py << "import numpy as np\n"
               << "from scipy.optimize import linprog, milp, LinearConstraint, Bounds\n"
               << "from scipy.sparse import lil_matrix\n"
               << "with open('" << inFile << "', 'r') as f:\n"
               << "    lines = [l.strip() for l in f if l.strip()]\n"
               << "p0 = lines[0].split()\n"
               << "nVars, nCons, sense, tlim = int(p0[0]), int(p0[1]), int(p0[2]), float(p0[3])\n"
               << "lb, ub, vtype, obj = np.zeros(nVars), np.zeros(nVars), np.zeros(nVars, dtype=int), np.zeros(nVars)\n"
               << "for i in range(nVars):\n"
               << "    a, b, t, c = lines[1+i].split()\n"
               << "    lb[i] = -np.inf if float(a) < -1e20 else float(a)\n"
               << "    ub[i] = np.inf if float(b) > 1e20 else float(b)\n"
               << "    vtype[i] = int(t)\n"
               << "    obj[i] = float(c) * sense\n"
               << "A = lil_matrix((nCons, nVars), dtype=float)\n"
               << "lhs, rhs = np.zeros(nCons), np.zeros(nCons)\n"
               << "for r in range(nCons):\n"
               << "    toks = lines[1+nVars+r].split()\n"
               << "    l_v, r_v, k = float(toks[0]), float(toks[1]), int(toks[2])\n"
               << "    lhs[r] = -np.inf if l_v < -1e20 else l_v\n"
               << "    rhs[r] = np.inf if r_v > 1e20 else r_v\n"
               << "    for j in range(k):\n"
               << "        A[r, int(toks[3+2*j])] = float(toks[3+2*j+1])\n"
               << "A_csc = A.tocsc()\n"
               << "is_mip = np.any(vtype > 0)\n"
               << "with open('" << outFile << "', 'w') as out:\n"
               << "    if not is_mip:\n"
               << "        eq = np.isclose(lhs, rhs)\n"
               << "        res = linprog(c=obj, A_ub=A_csc[~eq,:], b_ub=rhs[~eq], A_eq=A_csc[eq,:], b_eq=rhs[eq], bounds=list(zip(lb, ub)), method='highs')\n"
               << "        if res.success:\n"
               << "            out.write(f'1 {res.fun * sense:.10f} 0.0\\n')\n"
               << "            out.write(' '.join(f'{v:.10f}' for v in res.x) + '\\n')\n"
               << "        else:\n"
               << "            out.write('3 0.0 0.0\\n')\n"
               << "    else:\n"
               << "        integ = np.where(vtype == 2, 1, 0)\n"
               << "        if not np.any(integ): integ = np.where(vtype > 0, 1, 0)\n"
               << "        res = milp(c=obj, integrality=integ, bounds=Bounds(lb, ub), constraints=LinearConstraint(A_csc, lhs, rhs), options={'time_limit': tlim})\n"
               << "        if res.x is not None and len(res.x) > 0:\n"
               << "            st = 1 if res.status == 0 else 2\n"
               << "            gap = float(getattr(res, 'mip_gap', 0.0) or 0.0)\n"
               << "            out.write(f'{st} {res.fun * sense:.10f} {gap:.10f}\\n')\n"
               << "            out.write(' '.join(f'{v:.10f}' for v in res.x) + '\\n')\n"
               << "        else:\n"
               << "            out.write('3 0.0 0.0\\n')\n";
        }

        int ret = std::system("C:\\Python313\\python.exe _cplex_bridge_solver.py");
        std::remove(inFile);
        std::remove(pyFile);
        if (ret != 0) {
            std::remove(outFile);
            status_ = IloAlgorithm::Infeasible;
            return false;
        }

        std::ifstream out(outFile);
        if (!out.is_open()) {
            status_ = IloAlgorithm::Infeasible;
            return false;
        }
        int st = 3;
        out >> st >> objVal_ >> gap_;
        if (st == 1 || st == 2) {
            status_ = (st == 1) ? IloAlgorithm::Optimal : IloAlgorithm::Feasible;
            sol_.assign(nVars, 0.0);
            for (int i = 0; i < nVars; i++) out >> sol_[i];
            out.close();
            std::remove(outFile);
            return true;
        }
        out.close();
        std::remove(outFile);
        status_ = IloAlgorithm::Infeasible;
        return false;
    }
};

#endif /* ILOCPLEX_WINDOWS_COMPAT_H_ */
