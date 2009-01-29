/* vim: set sw=4 sts=4 et nofoldenable : */
#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>

//#include "NanoTimer.h"


class GLWidget : public QGLWidget
{
    Q_OBJECT

    public:
        GLWidget(QWidget * father = 0);
        ~GLWidget();

        QSize minimumSizeHint() const;
        QSize sizeHint() const;

        /*public slots:
          void setXRotation(int angle);
          void setYRotation(int angle);
          void setFlakeRec(int rec);
          void setBackgroundColor(int v);

signals:
void xRotationChanged(int angle);
void yRotationChanged(int angle);
*/
    protected:
        virtual void initializeGL();
        virtual void paintGL();
        virtual void resizeGL( int width, int height );
/*        virtual void mousePressEvent( QMouseEvent *event );
        virtual void mouseMoveEvent( QMouseEvent *event );
        virtual void timerEvent( QTimerEvent *event );
*/
    private:
/*        void makeSphere( unsigned int lati = 10, unsigned int longi = 10,
                const float radius = 1.0 ) const;
        void renderSphereFlake( const float scaling, unsigned int n_recursions ) const;

        void normalizeAngle(int *angle) const;*/
        void HsvToRgb( float h, float s, float v,
                float *r, float *g, float *b ) const;

        //void display_info( float framerate );

        GLuint m_object;
        int m_xRot, m_yRot, m_zRot;
        float m_xTrans, m_yTrans, m_zTrans;
        QPoint m_lastPos;
        float m_backgroundcolor[3];
        float m_curr_rot;
        unsigned int m_numFlakeRec;

        //NanoTimer m_timer;
};

#endif
