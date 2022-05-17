/**
******************************************************************************
* @file    main.c
* @author  SRA
* @version v2.5.0
* @date    23-Apr-2021
* @brief   Main program body
******************************************************************************
* @attention
*
* <h2><center>&copy; COPYRIGHT(c) 2016 STMicroelectronics</center></h2>
*
* Redistribution and use in source and binary forms, with or without modification,
* are permitted provided that the following conditions are met:
*   1. Redistributions of source code must retain the above copyright notice,
*      this list of conditions and the following disclaimer.
*   2. Redistributions in binary form must reproduce the above copyright notice,
*      this list of conditions and the following disclaimer in the documentation
*      and/or other materials provided with the distribution.
*   3. Neither the name of STMicroelectronics nor the names of its contributors
*      may be used to endorse or promote products derived from this software
*      without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
******************************************************************************
*/

/* Includes ------------------------------------------------------------------*/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "main.h"

#include "datalog_application.h"

/* FatFs includes component */
#include "ff_gen_drv.h"
#include "sd_diskio.h"

/* Private typedef -----------------------------------------------------------*/

/* Private define ------------------------------------------------------------*/

/* Private macro -------------------------------------------------------------*/
/* Private variables ---------------------------------------------------------*/
static volatile uint8_t MEMSInterrupt = 0;
static volatile uint8_t acquire_data_enable_request  = 1;
static volatile uint8_t acquire_data_disable_request = 0;

BSP_AUDIO_Init_t MicParams;

/* Private function prototypes -----------------------------------------------*/
void SW_Task1_Start(void);
void AudioProcess(void);
static uint8_t LSM6DSM_Sensor_IO_SPI_CS_Init(void);
static int32_t LSM6DSM_Sensor_IO_ITConfig( void );

/* Private functions ---------------------------------------------------------*/
volatile  int  index_buff=0;
volatile uint8_t write_flag=0;
volatile uint32_t MsInBuff=0;

uint16_t Audio_OUT_Buff[SIZE_BUFF];
volatile uint16_t cont=0;
uint16_t PCM_Buffer[AUDIO_CHANNELS * AUDIO_SAMPLING_FREQUENCY / 1000];



/**
* @brief  Main program
* @param  None
* @retval None
*/
int main( void )
{
  uint32_t msTick, msTickPrev = 0;  
  
  /* STM32L4xx HAL library initialization:
  - Configure the Flash prefetch, instruction and Data caches
  - Configure the Systick to generate an interrupt each 1 msec
  - Set NVIC Group Priority to 4
  - Global MSP (MCU Support Package) initialization
  */
  HAL_Init();
  
  /* Configure the system clock */
  SystemClock_Config();
  DATALOG_SD_Init();
  
  HAL_Delay(200);
  
  /* Initialize and disable all Chip Select */
  Sensor_IO_SPI_CS_Init_All();
  
  /* Configure and disable LSM6DSM Chip Select pins */
  LSM6DSM_Sensor_IO_SPI_CS_Init();
  
  /* Initialize and Enable the available sensors */
  MX_X_CUBE_MEMS1_Init();
  
  BSP_LED_Init(LED1);
  
  /* Configure Audio Input peripheral - DFSDM */  
  MicParams.BitsPerSample = 16;
  MicParams.ChannelsNbr = AUDIO_CHANNELS;
  MicParams.Device = AUDIO_IN_DIGITAL_MIC;
  MicParams.SampleRate = AUDIO_SAMPLING_FREQUENCY;
  MicParams.Volume = AUDIO_VOLUME_INPUT;
  
  BSP_AUDIO_IN_Init(BSP_AUDIO_IN_INSTANCE, &MicParams); 
  
  /* Start Microphone acquisition */
  BSP_AUDIO_IN_Record(BSP_AUDIO_IN_INSTANCE, (uint8_t *) PCM_Buffer, AUDIO_IN_BUFFER_SIZE);  
  
  while (1)
  {
    /* Check LSM6DSM Double Tap Event  */
    if(MEMSInterrupt)
    {
      MEMSInterrupt = 0;
      
      if(DoubleTap())
      {
        if (SD_Log_Enabled) 
        {
          SD_Log_Enabled=0;
          DATALOG_SD_Log_Disable();
          MsInBuff = 0;
        }
        else
        {
          while(SD_Log_Enabled != 1)
          {
            if(DATALOG_SD_Log_Enable())
            {
              SD_Log_Enabled=1;
            }
            else
            {
              DATALOG_SD_Log_Disable();
              DATALOG_SD_DeInit();
              DATALOG_SD_Init();
            }
          }
        }
      }
    }
    
    if (SD_Log_Enabled) 
    {
      msTick = HAL_GetTick();   /*led toggle*/
      if(msTick % 100 == 0 && msTickPrev != msTick)
      {
        msTickPrev = msTick;
        
        BSP_LED_Toggle(LED1);
      }
    }
    
    if(MsInBuff>0)
    {
      write_ms_on_sd();
      __disable_irq();
      MsInBuff--;
      __enable_irq();
    }
    
    /* Go to Sleep */
    if(MsInBuff==0)
    {
      __WFI();
    }
  }
}



void AudioProcess(void)
{
  /*for L4 PDM to PCM conversion is performed in hardware by DFSDM peripheral*/
  uint16_t index = 0;
  static uint16_t OUT_Buff_lvl = 0;
  
  if (SD_Log_Enabled)
  {
    for (index = 0; index < AUDIO_SAMPLING_FREQUENCY/1000 ; index++)
    {
      Audio_OUT_Buff[OUT_Buff_lvl] = PCM_Buffer[index];
      OUT_Buff_lvl = (OUT_Buff_lvl + 1)%SIZE_BUFF;
    }
    MsInBuff++;
  }
}


/**
* @brief  Half Transfer user callback, called by BSP functions.
* @param  None
* @retval None
*/
void BSP_AUDIO_IN_HalfTransfer_CallBack(uint32_t Instance)
{
  AudioProcess();
}

/**
* @brief  Transfer Complete user callback, called by BSP functions.
* @param  None
* @retval None
*/
void BSP_AUDIO_IN_TransferComplete_CallBack(uint32_t Instance)
{
  AudioProcess();
}


uint8_t LSM6DSM_Sensor_IO_SPI_CS_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct;
  
  /* Set all the pins before init to avoid glitch */
  HAL_GPIO_WritePin(BSP_LSM6DSM_CS_PORT, BSP_LSM6DSM_CS_PIN, GPIO_PIN_SET);
  
  GPIO_InitStruct.Speed = GPIO_SPEED_HIGH;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  
  BSP_LSM6DSM_CS_GPIO_CLK_ENABLE();
  GPIO_InitStruct.Pin = BSP_LSM6DSM_CS_PIN;
  HAL_GPIO_Init(BSP_LSM6DSM_CS_PORT, &GPIO_InitStruct);
  HAL_GPIO_WritePin(BSP_LSM6DSM_CS_PORT, BSP_LSM6DSM_CS_PIN, GPIO_PIN_SET);
  
  LSM6DSM_Sensor_IO_ITConfig();

  return BSP_ERROR_NONE;
}


/**
 * @brief  Configures sensor interrupts interface for LSM6DSM sensor.
 * @param  None
 * @retval BSP_ERROR_NONE in case of success
 * @retval COMPONENT_ERROR in case of failure
 */
int32_t LSM6DSM_Sensor_IO_ITConfig( void )
{
  /* At the moment this feature is only implemented for LSM6DSM */
  GPIO_InitTypeDef GPIO_InitStructureInt2;
  
  /* Enable INT2 GPIO clock */
  BSP_LSM6DSM_INT2_GPIO_CLK_ENABLE();
  
  /* Configure GPIO PINs to detect Interrupts */
  GPIO_InitStructureInt2.Pin = BSP_LSM6DSM_INT2;
  GPIO_InitStructureInt2.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStructureInt2.Speed = GPIO_SPEED_FAST;
  GPIO_InitStructureInt2.Pull  = GPIO_NOPULL;
  HAL_GPIO_Init(BSP_LSM6DSM_INT2_GPIO_PORT, &GPIO_InitStructureInt2);
  
  /* Enable and set EXTI Interrupt priority */
  HAL_NVIC_SetPriority(BSP_LSM6DSM_INT2_EXTI_IRQn, 0x08, 0x00);
  HAL_NVIC_EnableIRQ(BSP_LSM6DSM_INT2_EXTI_IRQn);
  
  return BSP_ERROR_NONE;
}


/**
* @brief  EXTI line detection callbacks
* @param  GPIO_Pin: Specifies the pins connected EXTI line
* @retval None
*/
void HAL_GPIO_EXTI_Callback( uint16_t GPIO_Pin )
{
  MEMSInterrupt=1;
}



#ifdef  USE_FULL_ASSERT

/**
* @brief  Reports the name of the source file and the source line number
*   where the assert_param error has occurred
* @param  file: pointer to the source file name
* @param  line: assert_param error line source number
* @retval None
*/
void assert_failed( uint8_t *file, uint32_t line )
{
  
  /* User can add his own implementation to report the file name and line number,
  ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  
  while (1)
  {}
}

#endif

/************************ (C) COPYRIGHT STMicroelectronics *****END OF FILE****/
