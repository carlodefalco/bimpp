#ifndef HAVE_MUPARSER_FUN_H
#define HAVE_MUPARSER_FUN_H

#include <muParser.h>
#include <memory>
#include <string>
#include <vector>
class muparser_fun
{
private :

  int N;
  std::vector<double> variables;
  mu::Parser p;

public :

  muparser_fun() :N(0) {};

  muparser_fun (const muparser_fun & m)
  {

    N = m.N;
    p = m.p;
    p.ClearVar();
    variables.clear();
    variables.resize(N);
    for(int i=0;i<N;i++)
    {
      std::string namevar = "x" + std::to_string(i);
      p.DefineVar (namevar, &variables[i]);
    }
  };


  muparser_fun& operator = (const muparser_fun &m)
  {
    N = m.N;
    p = m.p;
    p.ClearVar();
    variables.clear();
    variables.resize(N);
    for(int i=0;i<N;i++)
    {
      std::string namevar = "x" + std::to_string(i);
      p.DefineVar (namevar, &variables[i]);
    }
    return *this;
  }

  muparser_fun (const std::string & s, int N_):
  N(N_),variables(N_)
  {
    try
      {
        for(int i=0;i<N;i++)
        {
          std::string namevar = "x" + std::to_string(i);
          p.DefineVar (namevar, &variables[i]);//&var1);
          p.SetExpr (s);
        }
      }
    catch (mu::Parser::exception_type &e)
      { std::cerr << e.GetMsg () << std::endl; }
  };

  double operator() (std::vector<double> v)
  {
    double z;
    for(int i = 0; i < N; i++)
    {
        variables[i] = v[i];
    }
    try
      { z = p.Eval (); }
    catch (mu::Parser::exception_type &e)
      { std::cerr << e.GetMsg () << std::endl; }
    return (z);
  };

   muparser_fun & operator*=(const muparser_fun &rhs)
   {
     std::string prodexpr = "("+p.GetExpr()+")*("+rhs.p.GetExpr()+")";

     N = std::max(N,rhs.N);

     //p.ClearFun();
     p.ClearVar();

     variables.clear();
     variables.resize(N);

     for(int i=0;i<N;i++)
     {
       std::string namevar = "x" + std::to_string(i);
       p.DefineVar (namevar, &variables[i]);//&var1);
     }

     p.SetExpr (prodexpr);
     return *this;
   };

   muparser_fun & operator*=(const double &rhs)
   {
     std::string sumexpr = "("+p.GetExpr()+")*"+std::to_string(rhs);

     p.SetExpr (sumexpr);
     return *this;
   };



    muparser_fun & operator+=(const muparser_fun &rhs)
    {
      std::string sumexpr = p.GetExpr()+"+"+rhs.p.GetExpr();

      N = std::max(N,rhs.N);

      //p.ClearFun();
      p.ClearVar();

      variables.clear();
      variables.resize(N);

      for(int i=0;i<N;i++)
      {
        std::string namevar = "x" + std::to_string(i);
        p.DefineVar (namevar, &variables[i]);//&var1);
      }

      p.SetExpr (sumexpr);
      return *this;
    };


    muparser_fun & operator+=(const double &rhs)
    {
      std::string sumexpr = p.GetExpr()+"+"+std::to_string(rhs);

    //  p.ClearFun();

      p.SetExpr (sumexpr);
      return *this;
    };

    muparser_fun & operator-=(const muparser_fun &rhs)
    {
      std::string sumexpr = p.GetExpr()+"-"+rhs.p.GetExpr();

      N = std::max(N,rhs.N);

    //  p.ClearFun();
      p.ClearVar();

      variables.clear();
      variables.resize(N);

      for(int i=0;i<N;i++)
      {
        std::string namevar = "x" + std::to_string(i);
        p.DefineVar (namevar, &variables[i]);//&var1);
      }

      p.SetExpr (sumexpr);
      return *this;
    };


    const muparser_fun operator+(const muparser_fun &other) const
    {
      std::string sumexpr = p.GetExpr()+"+"+other.p.GetExpr();
      int Nnew = std::max(N,other.N);
      muparser_fun result(sumexpr, Nnew);
      return result;
    }

    const muparser_fun operator+(const double &other) const
    {
      std::string sumexpr = p.GetExpr()+"+"+std::to_string(other);
      muparser_fun result(sumexpr, N);
      return result;
    }

    const muparser_fun operator-(const muparser_fun &other) const
    {
      std::string sumexpr = p.GetExpr()+"-"+other.p.GetExpr();
      int Nnew = std::max(N,other.N);
      muparser_fun result(sumexpr, Nnew);
      return result;
    }

    const muparser_fun operator*(const muparser_fun &other) const
    {
      std::string prodexpr = "("+p.GetExpr()+")*("+other.p.GetExpr()+")";
      int Nnew = std::max(N,other.N);
      muparser_fun result(prodexpr, Nnew);
      return result;
    }

    const muparser_fun operator*(const double &other) const
    {
      std::string prodexpr = "("+p.GetExpr()+")*"+std::to_string(other);
      muparser_fun result(prodexpr, N);
      return result;
    }

    const muparser_fun operator/(const muparser_fun &other) const
        {
          std::string prodexpr = "("+p.GetExpr()+")/("+other.p.GetExpr()+")";
          int Nnew = std::max(N,other.N);
          muparser_fun result(prodexpr, Nnew);
          return result;
        }

        const muparser_fun operator/(const double &other) const
        {
          std::string prodexpr = "("+p.GetExpr()+")/"+std::to_string(other);
          muparser_fun result(prodexpr, N);
          return result;
        }

        std::vector<int> get_var_index()
        {
          auto var = p.GetUsedVar();
          std::vector<int> result;
          for(auto it = var.begin();it != var.end();it++)
          {
            std::string s = it->first;
            s = s.substr(1);
            result.push_back(std::stoi(s));
          }
          return result;
        }

        std::vector<int> get_var_index(int begin, int end)
        {
          auto var = p.GetUsedVar();
          std::vector<int> result;
          for(auto it = var.begin();it != var.end();it++)
          {
            std::string s = it->first;
            s = s.substr(1);
            int temp = std::stoi(s);
            if(temp >=begin && temp < end)
              result.push_back(temp);
          }
          return result;
        }

        muparser_fun & operator/=(const muparser_fun &rhs)
   {
     std::string prodexpr = "("+p.GetExpr()+")/("+rhs.p.GetExpr()+")";

     N = std::max(N,rhs.N);

     //p.ClearFun();
     p.ClearVar();

     variables.clear();
     variables.resize(N);

     for(int i=0;i<N;i++)
     {
       std::string namevar = "x" + std::to_string(i);
       p.DefineVar (namevar, &variables[i]);//&var1);
     }

     p.SetExpr (prodexpr);
     return *this;
   };

   muparser_fun & operator/=(const double &rhs)
   {
     std::string sumexpr = "("+p.GetExpr()+")*"+std::to_string(rhs);

     p.SetExpr (sumexpr);
     return *this;
   };




};

#endif
