#ifndef AUXGRID
#define AUXGRID
#include "Include.hh"

enum BoundaryType{Neumann, Dirichlet, Periodic, Undefined};

struct Boundary{
	BoundaryType type;
	Grid value;
    bool ghost;

	Boundary(){ type = Undefined; ghost = false;}

	Boundary& operator<< (const BoundaryType &_type){ type = _type; return *this;}
	Boundary& operator| (const Grid &_value){ value = _value; return *this;}

	friend std::ostream& operator<< (std::ostream &os, const Boundary &other){	
		os << "\tType : ";
		switch (other.type){
			case Neumann:
				os << "Neumann  \t";
				break;
			case Dirichlet:
				os << "Dirichlet\t";
				break;
			case Periodic:
				os << "Periodic \t";
				break;
			case Undefined:
				os << "Undefined\t";
				break;
		}
		os << " | ";
		if (other.value.size() > 4){
			os << "[" << other.value(0) << ", "
				<< other.value(1) << ", "
				<< "... "
				<< other.value(other.value.size()-2) << ", "
				<< other.value(other.value.size()-1) << "] ";
		} else{
            if (other.value.rows() == 1) os << other.value;
            else os << other.value.transpose();
		}
		return os;
	}
};

class AuxGrid{
public:
    Grid left, right, bottom, top;

	friend std::ostream& operator<< (std::ostream &os, const AuxGrid &other){	
		os << "Left" << other.l << std::endl;
		os << "Right"<< other.r << std::endl;
		os << "Bottom"<< other.b << std::endl;
		os << "Top" << other.t;
		return os;
	}
    void set_row(BoundaryType type, Grid value){r = l << type | value;}
    void set_row(BoundaryType type){r = l << type;}
    void set_col(BoundaryType type, Grid value){b = t << type | value;}
    void set_col(BoundaryType type){b = t << type;}
    void set_left(BoundaryType type, Grid value){ l << type | value;}
    void set_left(BoundaryType type){ l << type;}
    void set_right(BoundaryType type, Grid value){ r << type | value;}
    void set_right(BoundaryType type){ r << type;}
    void set_bottom(BoundaryType type, Grid value){ b << type | value;}
    void set_bottom(BoundaryType type){ b << type;}
    void set_top(BoundaryType type, Grid value){ t << type | value;}
    void set_top(BoundaryType type){ t << type;}
    void set_all(BoundaryType type, Grid value){r = l = b = t << type | value;}
    void set_all(BoundaryType type){r = l = b = t << type;}

    void set_row_ghost() {l.ghost = true; r.ghost = true;}
    void set_col_ghost() {b.ghost = true; t.ghost = true;}
    bool row_ghost() {return (l.ghost && r.ghost);}
    bool col_ghost() {return (b.ghost && t.ghost);}
    bool row_periodic() {return (l.type == Periodic && r.type == Periodic);}
    bool col_periodic() {return (b.type == Periodic && t.type == Periodic);}
    void update(const Grid &var){
        switch (l.type){
            case Neumann:
                if (l.value.size() == 0) l.value = ZERO2(1, var.cols());
                left = var.row(0) - l.value;
                break;
            case Dirichlet:
                if (l.value.size() == 0) l.value = var.row(0);
                left = l.value;
                break;
            case Periodic:
                left = var.row(var.rows() - 1);
                break;
            case Undefined:
                left = ZERO2(1, var.cols());
                break;
        }
        switch (r.type){
            case Neumann:
                if (r.value.size() == 0) r.value = ZERO2(1, var.cols());
                right = var.row(var.rows() - 1) - r.value;
                break;
            case Dirichlet:
                if (r.value.size() == 0) r.value = var.row(var.rows() - 1);
                right = r.value;
                break;
            case Periodic:
                right = var.row(0);
                break;
            case Undefined:
                right = ZERO2(1, var.cols());
                break;
        }
        switch (b.type){
            case Neumann:
                if (b.value.size() == 0) b.value = ZERO2(var.rows(), 1);
                bottom = var.col(0) - b.value;
                break;
            case Dirichlet:
                if (b.value.size() == 0) b.value = var.col(0);
                bottom = b.value;
                break;
            case Periodic:
                bottom = var.col(var.cols() - 1);
                break;
            case Undefined:
                bottom = ZERO2(var.rows(), 1);
                break;
        }
        switch (t.type){
            case Neumann:
                if (t.value.size() == 0) t.value = ZERO2(var.rows(), 1);
                top = var.col(var.cols() - 1) - t.value;
                break;
            case Dirichlet:
                if (t.value.size() == 0) t.value = var.col(var.cols() - 1);
                top = t.value;
                break;
            case Periodic:
                top = var.col(0);
                break;
            case Undefined:
                top = ZERO2(var.rows(), 1);
                break;
        }
    }
    AuxGrid operator*(const AuxGrid &aux){
        AuxGrid result;
        result.left = this->left * aux.left;
        result.right= this->right * aux.right;
        result.bottom= this->bottom * aux.bottom;
        result.top= this->top * aux.top;
        return result;
    }
protected:
    Boundary l, r, b, t;
};
#endif
