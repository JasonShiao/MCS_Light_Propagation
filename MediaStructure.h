#ifndef MEDIASTRUCTURE_H_INCLUDED
#define MEDIASTRUCTURE_H_INCLUDED

class Media{
   public:
      Media();

      double get_mu_a( double z){ return mu_a; };
      double get_length(){ return length; };

      virtual ~Media();
   protected:
   private:
      double mu_a = 10;
      double length = 1;

};

#endif // MEDIASTRUCTURE_H_INCLUDED
